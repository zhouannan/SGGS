
#pragma once

#include <string>
#include <vector>

#include "Vector3.h"
#include "Attributes.h"
#include "Monitor.h"

using std::string;
using std::vector;

class Source;
class State;

namespace chunker_countsort_laszip
{

	struct ClusterInfo
	{
		Vector3 min;
		Vector3 max;
		int count;
		double dis;
	};

	void doChunking(vector<Source> sources, string targetDir, Vector3 min, Vector3 max, State &state, 
	Attributes outputAttributes, Monitor *monitor, Options options,std::unordered_map<int, chunker_countsort_laszip::ClusterInfo> clusterInfoMap);

}