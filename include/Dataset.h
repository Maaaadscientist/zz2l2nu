#ifndef DATASET_H_
#define DATASET_H_

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <yaml-cpp/yaml.h>

#include <TChain.h>
#include <TTreeReader.h>

#include <Options.h>


/**
 * \brief Input files included in a dataset and metadata about it
 * 
 * An object of this class is created from a text file definining a dataset.
 * Several file formats are supported, as described
 * <a href="https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/wikis/dataset-definitions">here</a>:
 * -# Full YAML dataset definition file.
 * -# YAML dataset definition fragment that provides several parameters and
 *    points to a stem definition fragment. The stem fragment is included into
 *    it. This class searches the stem fragment within YAML files whose
 *    locations are specified in the master configuration.
 *
 * The access to the list of input files and some generic parameters is
 * provided with dedicated methods. All other parameters are accessible via
 * YAML::Node returned by \ref Parameters.
 */
class DatasetInfo {
 public:
  /// Constructor from a dataset definition file
  DatasetInfo(std::filesystem::path const &path, Options const &options);

  /// Returns path to dataset definition file used
  std::filesystem::path const &DefinitionFile() const {
    return definitionFile_;
  }

  /**
   * \brief Returns the cross section, in pb
   *
   * Undefined for real data.
   */
  double CrossSection() const {
    return crossSection_;
  }

  /**
   * \brief Returns paths to all input files in the dataset
   *
   * The paths can include protocol prefix for ROOT::TFile::Open.
   */
  std::vector<std::string> const &Files() const {
    return files_;
  }

  /// Indicates whether this is simulation or real data
  bool IsSimulation() const {
    return isSimulation_;
  };

  /**
   * \brief Returns mean nominal weight in the dataset
   *
   * Undefined for real data.
   */
  double MeanWeight() const {
    return meanWeight_;
  }
  
  /// Returns unique name of the dataset
  std::string const &Name() const {
    return name_;
  }

  /**
   * \brief Returns total number of events before any selection
   *
   * Undefined for real data.
   */
  int64_t NumEventsTotal() const {
    return numEventsTotal_;
  }

  /**
   * \brief Parameters extracted from dataset definition file
   *
   * This includes all parameters except the paths to input files.
   */
  YAML::Node const &Parameters() const {
    return parameters_;
  }

 private:
  /// Finds stem dataset definition with given name
  YAML::Node const FindStem(std::string_view name) const;

  /// Get a scalar node from dataset definition file with error reporting
  YAML::Node const GetNode(YAML::Node const root, std::string const &key) const;

  /**
   * \brief Reads YAML dataset definition file
   *
   * The stem definition is incorporated if needed.
   */
  void ReadYaml(std::filesystem::path const &path);

  /**
   * \brief Incorporate stem definition into the given dataset definition
   * fragment
   *
   * The given node is modified in place. Key "stem" is removed.
   */
  void SpliceYaml(YAML::Node info) const;

  /**
   * \brief Locations of files with dataset stems
   *
   * Only needed to construct the full dataset definition from a YAML fragment.
   */
  std::vector<std::filesystem::path> stemsFiles_;

  /**
   * \brief Path to dataset definition file
   *
   * Mostly needed for error reports.
   */
  std::filesystem::path definitionFile_;

  /**
   * \brief Paths to input files in the dataset
   *
   * These can include protocol prefix for ROOT::TFile::Open and therefore
   * cannot be represented with std::filesystem::path.
   */
  std::vector<std::string> files_;

  /**
   * \brief Parameters extracted from dataset definition file
   *
   * See \ref Parameters.
   */
  YAML::Node parameters_;

  /// Name that uniquely identifies the dataset
  std::string name_;

  /// Indicates whether this is simulation or real data
  bool isSimulation_;

  /**
   * \brief Cross section, in pb
   *
   * Undefined for real data.
   */
  double crossSection_;

  /**
   * \brief Total number of events before any selection
   *
   * Undefined for real data.
   */
  int64_t numEventsTotal_;

  /**
   * \brief Mean nominal event weight
   *
   * Undefined for real data.
   */
  double meanWeight_;
};


/**
 * \brief Interface to read a dataset
 *
 * This class is an aggregation of a DatasetInfo and a TTreeReader. The latter
 * one implements reading of input files included in the dataset. Reading of a
 * subset of files in the dataset is supported, as needed for parallel
 * processing.
 *
 * When a consumer registers a new branch to read, this modifies the Dataset
 * object (meaning that \ref Reader, which has to be called for this, is not a
 * constant method). This can be thought of as a (dynamic) change in the content
 * of the dataset.
 */
class Dataset {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] info       DatasetInfo object that describes the dataset.
   * \param[in] skipFiles  Number of input files in the dataset to skip,
   *   starting from the beginning.
   * \param[in] maxFiles   Number of input files to include, starting after the
   *   skipped ones. Special value -1 indicates that all remaining files need to
   *   be included.
   */
  Dataset(DatasetInfo info, int skipFiles = 0, int maxFiles = -1);

  /// Returns associated DatasetInfo object
  DatasetInfo const &Info() const {
    return info_;
  }

  /**
   * \brief Sets the next entry in the dataset as the current one
   *
   * \return False if the previous entry was the last one in the dataset; true
   * otherwise.
   */
  bool NextEntry() {
    return reader_.Next();
  }

  /// Returns number of entries in the selected input files from the dataset
  int64_t NumEntries() {
    return reader_.GetEntries(true);
  }

  /// Returns reader for selected files
  TTreeReader &Reader() {
    return reader_;
  }

  /// Returns paths to selected input files
  std::vector<std::string> const &SelectedFiles() const {
    return selectedFiles_;
  }

  /// Sets the current entry for the reader
  void SetEntry(int64_t index) {
    reader_.SetEntry(index);
  }

 private:
  /// Associated DatasetInfo object
  DatasetInfo info_;

  /// Selected files from the dataset
  std::vector<std::string> selectedFiles_;

  /// TChain containing selected files
  TChain chain_;

  /// Reader associated with \ref chain_
  TTreeReader reader_;
};

#endif  // DATASET_H_

