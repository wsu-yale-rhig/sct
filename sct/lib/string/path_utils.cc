#include "sct/lib/string/path_utils.h"

namespace sct {

    string GetFileName(const string& path) {
        const char separator = '/';
        size_t position = path.rfind(separator);
        if (position != string::npos)
            return path.substr(position + 1, string::npos);
        else
            return path;
    }


    string GetPath(const string& path) {
        const char separator = '/';
        size_t position = path.rfind(separator);
        if (position != string::npos)
            return path.substr(0, position);
        else 
            return path;
    }

} // namespace sct
