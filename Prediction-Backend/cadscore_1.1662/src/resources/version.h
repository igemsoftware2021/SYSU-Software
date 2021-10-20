#ifndef RESOURCES_VERSION_H_
#define RESOURCES_VERSION_H_

#include <string>

namespace resources
{

inline std::string get_version_string()
{
	static std::string version_string("cadscore_1.1662");
	return version_string;
}

}

#endif /* RESOURCES_VERSION_H_ */
