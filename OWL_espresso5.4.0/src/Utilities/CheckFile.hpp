#include <sys/stat.h>
#include <string>
#include <fstream>

inline bool file_exists (const std::string& filename) {
  std::ifstream f(filename);
  //std::ifstream f(filename.c_str());
  return f.good();
}

inline bool file_exists (const char* filename) {
  struct stat buffer;  
   return (stat (filename, &buffer) == 0);
}

