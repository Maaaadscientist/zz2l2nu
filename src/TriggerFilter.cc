#include <TriggerFilter.h>

#include <limits>

#include <Logger.h>
#include <HZZException.h>


TriggerFilter::TriggerFilter(
    Dataset &dataset, Options const &options, RunSampler const *runSampler)
    : runSampler_{runSampler}, cache_{dataset.Reader()},
      reader_{&dataset.Reader()}, treeIndex_{-1} {
  if (not reader_->IsChain())
    throw HZZException(
        "TriggerFilter expects that the dataset contains a TChain.");
  chain_ = dynamic_cast<TChain *>(reader_->GetTree());

  auto const config = options.GetConfig()["trigger_filter"];
  if (not config)
      throw HZZException(
          "Section \"trigger_filter\", which is required by TriggerFilter, is "
          "missing in the master configuration.");
  LoadConfig(config);

  for (auto const &[channelName, channel] : channels_) {
    LOG_TRACE << "Triggers in channel \"" << channelName << "\"";
    for (auto const &trigger : channel.triggers)
      LOG_TRACE << "  " << trigger.trigger->name << " [" << trigger.minRun
          << ", " << trigger.maxRun << "]\n";
  }
}


bool TriggerFilter::operator()(std::string_view channel) const {
  if (cache_.IsUpdated())
    Build();
  auto const res = channels_.find(channel);
  if (res == channels_.end()) {
    HZZException exception;
    exception << "Unknown trigger channel \"" << channel << "\" requested.";
    throw exception;
  }
  return res->second.decision;
}


bool TriggerFilter::TriggerInPeriod::GetDecision(run_t run) const {
  if (run < minRun or run > maxRun)
    return false;
  return (trigger->branch and trigger->decision);
}


bool TriggerFilter::Channel::Collect(run_t run) const {
  decision = false;
  for (auto const trigger : triggers) {
    if (trigger.GetDecision(run)) {
      decision = true;
      break;
    }
  }
  return decision;
}


void TriggerFilter::Build() const {
  // Update statuses and addresses of trigger branches if the input chain has
  // switched to a new tree
  int const curTreeIndex = chain_->GetTreeNumber();
  if (curTreeIndex != treeIndex_) {
    treeIndex_ = curTreeIndex;
    for (auto const &trigger : triggers_) {
      trigger.branch = chain_->GetBranch(("HLT_" + trigger.name).c_str());
      if (trigger.branch)
        trigger.branch->SetAddress(&trigger.decision);
    }
  }

  // Read trigger decisions into the buffers
  int64_t const entryCurTree = reader_->GetCurrentEntry()
      - chain_->GetTreeOffset()[treeIndex_];
  for (auto const &trigger : triggers_) {
    if (trigger.branch)
      trigger.branch->GetEntry(entryCurTree);
  }

  // Collect trigger decisions for all channels
  for (auto const &[name, channel] : channels_)
    channel.Collect(runSampler_->Get());
}


void TriggerFilter::LoadConfig(YAML::Node const config) {
  if (not config.IsMap())
    throw HZZException("Node [trigger_filter] must be a mapping.");

  for (auto const &channelNode : config) {
    Channel channel;
    auto const channelName = channelNode.first.as<std::string>();
    auto const &channelConfig = channelNode.second;
    if (not channelConfig.IsSequence() or channelConfig.size() == 0) {
      HZZException exception;
      exception << "Node [trigger_filter][" << channelName
          << "] must be a non-empty sequence.";
      throw exception;
    }

    for (int blockIndex = 0; blockIndex < int(channelConfig.size());
         ++blockIndex) {
      auto const &triggerBlockNode = channelConfig[blockIndex];
      auto const &triggerNamesNode = triggerBlockNode["triggers"];
      if (not triggerNamesNode) {
        HZZException exception;
        exception << "Mandatory node [triggers] is missing in [trigger_filter]["
            << channelName << "][" << blockIndex << "].";
        throw exception;
      }
      if (not triggerNamesNode.IsSequence() or triggerNamesNode.size() == 0) {
        HZZException exception;
        exception << "Node [trigger_filter][" << channelName << "]["
            << blockIndex << "][triggers] must be a non-empy sequence";
        throw exception;
      }
      auto const triggerNames = triggerNamesNode.as<std::vector<std::string>>();

      run_t minRun, maxRun;
      auto const &runRangeNode = triggerBlockNode["run_range"];
      if (not runRangeNode) {
        minRun = std::numeric_limits<run_t>::min();
        maxRun = std::numeric_limits<run_t>::max();
      } else {
        if (not runRangeNode.IsSequence() or runRangeNode.size() != 2) {
          HZZException exception;
          exception << "Node [trigger_filter][" << channelName << "]["
              << blockIndex << "][run_range] must be a sequence of two "
              << "elements.";
          throw exception;
        }
        minRun = runRangeNode[0].as<run_t>();
        maxRun = runRangeNode[1].as<run_t>();
        if (minRun > maxRun) {
          HZZException exception;
          exception << "Incorrect ordering in [trigger_filter][" << channelName
              << "][" << blockIndex << "][run_range].";
          throw exception;
        }
      }

      for (auto const &triggerName : triggerNames) {
        auto const res = triggers_.emplace(triggerName);
        Trigger const *trigger = &*res.first;
        channel.triggers.emplace_back(trigger, minRun, maxRun);
      }
    }

    channels_[channelName] = channel;
  }
}
