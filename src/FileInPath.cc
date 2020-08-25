#include <FileInPath.h>

#include <HZZException.h>

#include <cstdlib>

#include <boost/algorithm/string/predicate.hpp>


namespace fs = std::filesystem;


void FileInPath::AddLocation(fs::path const &path) {
  GetInstance().locations_.emplace_back(path);
}


fs::path FileInPath::Resolve(fs::path const &subDir, fs::path const &path) {

  // Check if an absolute or explicit relative path is provided. In this case
  // the first argument is ignored and the path is returned as is.
  auto const pathString = path.string();

  if (boost::starts_with(pathString, "/") or
      boost::starts_with(pathString, "./") or
      boost::starts_with(pathString, "../")) {

    // Make sure the requested file exists
    if (not fs::exists(path)) {
      HZZException exception;
      exception << "FileInPath::Resolve: "
        "Requested file with absolute or explicit relative path " <<
        path << " does not exist or is not a regular file.";
      throw exception;
    }

    return path;
  }


  // Loop over all possible locations, giving preference to ones added later
  auto const &locations = GetInstance().locations_;

  for (auto locationIt = locations.crbegin(); locationIt != locations.crend();
      ++locationIt) {
    fs::path tryPath;

    // Try to resolve the path using the provided subdirectory
    if (not subDir.empty()) {
      tryPath = *locationIt / subDir / path;

      if (fs::exists(tryPath))
        return tryPath;
    }

    // Try to resolve the path ignoring the subdir
    tryPath = *locationIt / path;

    if (fs::exists(tryPath))
      return tryPath;
  }


  // If none of the above succeeds, try to resolve the path with respect to the
  // currect working directory, with and without the subdirectory
  if (not subDir.empty()) {
    fs::path tryPath = subDir / path;

    if (fs::exists(tryPath))
      return fs::current_path() / tryPath;
  }

  if (fs::exists(path))
    return fs::current_path() / path;


  // If the workflow has reached this point, the path has not been resolved
  HZZException exception;
  exception << "FileInPath::Resolve: Failed to resolve path " << path <<
    " within (optional) subdirectory " << subDir <<
    ". Following locations have been tried (with and without the "
    "subdirectory):\n";

  for (auto locationIt = locations.crbegin(); locationIt != locations.crend();
       ++locationIt)
    exception << "  " << *locationIt << '\n';

  exception << "  " << fs::current_path() << '\n';
  throw exception;
}


fs::path FileInPath::Resolve(fs::path const &path) {
  return Resolve({}, path);
}


FileInPath::FileInPath() {
  // Read the install path from the environment
  char const *installPath = std::getenv("HZZ2L2NU_BASE");

  if (not installPath)
    throw HZZException(
      "FileInPath::FileInPath: Mandatory environment variable "
      "HZZ2L2NU_BASE is not defined.");


  // Specify default locations. Cannot use method AddLocation here since it will
  // indirectly call the constructor.
  locations_.emplace_back(fs::path{installPath} / "config");
  locations_.emplace_back(fs::path{installPath} / "data");
}


FileInPath &FileInPath::GetInstance() {
  static FileInPath instance;
  return instance;
}
