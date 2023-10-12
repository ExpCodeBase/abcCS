#include <time.h>  /* gettimeofday */
#include "utility.h"

/* Returns current time in seconds as a double value
and the precision can be as accurate as microseconds
(10^-6 second)
*/
double get_cur_time() {
	return 0;
}


std::ostream& operator<<(std::ostream& out, const Index_update value) {
	switch (value)
	{
		case Index_update::withlimit: return out << "withlimit";
		case Index_update::withlimit_base_opt: return out << "withlimit_base_opt";
		case Index_update::withlimit_dfs: return out << "withlimit_dfs";
		case Index_update::withlimit_dfs_parallel: return out << "withlimit_dfs_parallel";
		case Index_update::withlimit_parallel: return out << "withlimit_parallel";
		case Index_update::withoutlimit: return out << "withoutlimit";
		case Index_update::batch: return out << "batch";
		case Index_update::batch_parallel: return out << "batch_parallel";
	};
	return out << static_cast<std::uint16_t>(value);
}

std::ostream& operator<<(std::ostream& out, const Datasets value) {
	switch (value)
	{
		case Datasets::RE: return out << "RE";
		case Datasets::TR: return out << "TR";
		case Datasets::DUI: return out << "DUI";
		
	};
	return out << static_cast<std::uint16_t>(value);
}

