#include "Logger.h"

#include <unistd.h>
#include <cstdio>
#include <iostream>

#include <boost/core/null_deleter.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>


// Format function for log record
void FormatRecord(boost::log::record_view const &record,
                  boost::log::formatting_ostream &stream) {
  auto const severity = boost::log::extract<Logger::SeverityLevel>(
    "Severity", record);

  if (severity) {
    switch (severity.get()) {
      case Logger::SeverityLevel::kWarning:
        stream << "[WARN] ";
        break;

      case Logger::SeverityLevel::kError:
        stream << "[ERROR] ";
        break;

      default:
        break;
    }
  }

  stream << record[boost::log::expressions::smessage];
}


// Format function for log record, which uses coloured output
void FormatRecordColour(boost::log::record_view const &record,
                        boost::log::formatting_ostream &stream) {
  auto const severity = boost::log::extract<Logger::SeverityLevel>(
    "Severity", record);
  bool needReset = false;

  if (severity) {
    switch (severity.get()) {
      case Logger::SeverityLevel::kWarning:
        stream << "\033[33m[WARN] ";  // Yellow
        needReset = true;
        break;

      case Logger::SeverityLevel::kError:
        stream << "\033[31m[ERROR] ";  // Red
        needReset = true;
        break;

      default:
        break;
    }
  }

  stream << record[boost::log::expressions::smessage];

  if (needReset)
    stream << "\033[0m";
}


Logger::Logger() {
  // Construct and register a sink that writes output to stderr
  boost::shared_ptr<std::ostream> clogStream(&std::clog, boost::null_deleter());
  sink_ = boost::make_shared<decltype(sink_)::element_type>();
  sink_->locked_backend()->add_stream(clogStream);
  boost::log::core::get()->add_sink(sink_);

  // Flush associated stream after each log record
  sink_->locked_backend()->auto_flush(true);

  // Set formatter. Use coloured output if running in an interactive terminal.
  if (isatty(fileno(stderr)))
    sink_->set_formatter(&FormatRecordColour);
  else
    sink_->set_formatter(&FormatRecord);

  // By default, use the lowest severity threshold for reporting
  auto severity = boost::log::expressions::attr<SeverityLevel>("Severity");
  sink_->set_filter(severity >= SeverityLevel::kTrace);
}


void Logger::SetLevel(SeverityLevel threshold) {
  auto severity = boost::log::expressions::attr<SeverityLevel>("Severity");
  GetInstance().sink_->set_filter(severity >= threshold);
}

