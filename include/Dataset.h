#ifndef DATASET_H_
#define DATASET_H_

#include <cstdint>
#include <filesystem>
#include <vector>

#include <TChain.h>
#include <TTreeReader.h>


/// Input files included in a dataset and metadata about it
class DatasetInfo {
 public:
  /// Constructor from a dataset definition file
  DatasetInfo(std::filesystem::path const &path);

  /// Returns path to dataset definition file used
  std::filesystem::path const &DefinitionFile() const {
    return definitionFile_;
  }

  /// Returns paths to all input files in the dataset
  std::vector<std::filesystem::path> const &Files() const {
    return files_;
  }

  /// Indicates whether this is simulation or real data
  bool IsSimulation() const {
    return isSimulation_;
  };

 private:
  /// Parses txt dataset definition file
  void ParseText(std::filesystem::path const &path);

  /**
   * \brief Path to dataset definition file
   *
   * Mostly needed for error reports.
   */
  std::filesystem::path definitionFile_;

  /// Paths to input files in the dataset
  std::vector<std::filesystem::path> files_;

  /// Indicates whether this is simulation or real data
  bool isSimulation_;
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
  int64_t NumEntries() const {
    return reader_.GetEntries(true);
  }

  /// Returns reader for selected files
  TTreeReader &Reader() {
    return reader_;
  }

  /// Returns paths to selected input files
  std::vector<std::filesystem::path> const &SelectedFiles() const {
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
  std::vector<std::filesystem::path> selectedFiles_;

  /// TChain containing selected files
  TChain chain_;

  /// Reader associated with \ref chain_
  TTreeReader reader_;
};

#endif  // DATASET_H_

