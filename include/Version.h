#ifndef VERSION_H_
#define VERSION_H_

#include <string>


/**
 * \brief Provides access to the version of the package
 *
 * Implemented as a Mayer's singleton.
 */
class Version {
 public:
  /// Returns the hash of the Git commit
  static std::string const &Commit() {
    return GetInstance().commit_;
  }

 private:
  Version();

  Version(Version const &) = delete;
  Version &operator=(Version const &) = delete;

  /// Returns the only instance of this class
  static Version &GetInstance() {
    static Version instance;
    return instance;
  }
  
  /// Hash of the commit
  std::string commit_;
};

#endif  // VERSION_H_

