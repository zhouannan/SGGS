
#include <iostream>
#include <filesystem>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <memory>
#include <atomic>

#include "chunker_countsort_laszip.h"

#include "Attributes.h"
#include "converter_utils.h"
#include "unsuck/unsuck.hpp"
#include "unsuck/TaskPool.hpp"
#include "Vector3.h"
#include "ConcurrentWriter.h"

#include "json/json.hpp"
#include "laszip/laszip_api.h"
#include "LasLoader/LasLoader.h"
#include "PotreeConverter.h"
#include "logger.h"

using json = nlohmann::json;

using std::atomic_int32_t;
using std::cout;
using std::endl;
using std::lock_guard;
using std::make_shared;
using std::mutex;
using std::shared_ptr;
using std::thread;
using std::unique_ptr;
using std::unordered_map;

namespace fs = std::filesystem;

namespace chunker_countsort_laszip
{

	auto numChunkerThreads = getCpuData().numProcessors;
	auto numFlushThreads = getCpuData().numProcessors;

	// auto numChunkerThreads = 1;
	// auto numFlushThreads = 1;

	int maxPointsPerChunk = 5000000;
	int gridSize = 128;
	mutex mtx_attributes;

	struct Point
	{
		double x;
		double y;
		double z;
	};

	struct ThreadSafeCluster
	{
		std::unordered_map<int, int> clusterMap;
		std::mutex mtx; // 互斥锁

		ThreadSafeCluster() = default;
	};

	struct CountResult
	{
		vector<std::atomic_int32_t> grid;
		vector<std::shared_ptr<ThreadSafeCluster>> gridCluster;
	};

	unordered_map<int, vector<vector<Point>>> buckets;

	ConcurrentWriter *writer = nullptr;

	struct Node
	{

		string id = "";
		int64_t level = 0;
		int64_t x = 0;
		int64_t y = 0;
		int64_t z = 0;
		int64_t size;
		int64_t numPoints;
		std::unordered_map<int, int> clusterMap;

		Node(string id, int numPoints)
		{
			this->id = id;
			this->numPoints = numPoints;
		}

		Node(string id, int numPoints, std::unordered_map<int, int> clusterMap)
		{
			this->id = id;
			this->numPoints = numPoints;
			this->clusterMap = clusterMap;
		}
	};

	vector<Node> nodes;

	string toNodeID(int level, int gridSize, int64_t x, int64_t y, int64_t z)
	{

		string id = "r";

		int currentGridSize = gridSize;
		int lx = x;
		int ly = y;
		int lz = z;

		for (int i = 0; i < level; i++)
		{

			int index = 0;

			if (lx >= currentGridSize / 2)
			{
				index = index + 0b100;
				lx = lx - currentGridSize / 2;
			}

			if (ly >= currentGridSize / 2)
			{
				index = index + 0b010;
				ly = ly - currentGridSize / 2;
			}

			if (lz >= currentGridSize / 2)
			{
				index = index + 0b001;
				lz = lz - currentGridSize / 2;
			}

			id = id + to_string(index);
			currentGridSize = currentGridSize / 2;
		}

		return id;
	}

	template <class T>
	double asDouble(uint8_t *data)
	{
		T value = reinterpret_cast<T *>(data)[0];

		return double(value);
	}

	void countGridPointNumber(double *coordinates, vector<std::atomic_int32_t> &grid, double dGridSize,
							  Vector3 posOffset, Vector3 posScale, Vector3 min, Vector3 max, Vector3 size, vector<std::shared_ptr<ThreadSafeCluster>> &gridCluster, int label)
	{
		double x = coordinates[0];
		double y = coordinates[1];
		double z = coordinates[2];

		float_t X = float_t((x - posOffset.x) / posScale.x);
		float_t Y = float_t((y - posOffset.y) / posScale.y);
		float_t Z = float_t((z - posOffset.z) / posScale.z);

		double ux = (double(X) * posScale.x + posOffset.x - min.x) / size.x;
		double uy = (double(Y) * posScale.y + posOffset.y - min.y) / size.y;
		double uz = (double(Z) * posScale.z + posOffset.z - min.z) / size.z;

		bool inBox = ux >= 0.0 && uy >= 0.0 && uz >= 0.0;
		inBox = inBox && ux <= 1.01 && uy <= 1.01 && uz <= 1.01;

		if (!inBox)
		{
			stringstream ss;
			ss << "encountered point outside bounding box." << endl;
			ss << "box.min: " << min.toString() << endl;
			ss << "box.max: " << max.toString() << endl;
			ss << "point: " << Vector3(x, y, z).toString() << endl;
			ss << "PotreeConverter requires a valid bounding box to operate." << endl;
			ss << "Please try to repair the bounding box, e.g. using lasinfo with the -repair_bb argument." << endl;
			logger::ERROR(ss.str());

			exit(123);
		}

		int64_t ix = int64_t(std::min(dGridSize * ux, dGridSize - 1.0));
		int64_t iy = int64_t(std::min(dGridSize * uy, dGridSize - 1.0));
		int64_t iz = int64_t(std::min(dGridSize * uz, dGridSize - 1.0));

		int64_t index = ix + iy * gridSize + iz * gridSize * gridSize;

		grid[index]++;
		std::shared_ptr<ThreadSafeCluster> &cluster = gridCluster[index];
		std::lock_guard<std::mutex> lock(cluster->mtx); // 锁住互斥锁
		cluster->clusterMap[label]++;					// 安全地修改 unordered_map
	};

	void countGridPointNumber(double *coordinates, vector<std::atomic_int32_t> &grid, double dGridSize,
							  Vector3 posOffset, Vector3 posScale, Vector3 min, Vector3 max, Vector3 size)
	{
		double x = coordinates[0];
		double y = coordinates[1];
		double z = coordinates[2];

		float_t X = float_t((x - posOffset.x) / posScale.x);
		float_t Y = float_t((y - posOffset.y) / posScale.y);
		float_t Z = float_t((z - posOffset.z) / posScale.z);

		double ux = (double(X) * posScale.x + posOffset.x - min.x) / size.x;
		double uy = (double(Y) * posScale.y + posOffset.y - min.y) / size.y;
		double uz = (double(Z) * posScale.z + posOffset.z - min.z) / size.z;

		bool inBox = ux >= 0.0 && uy >= 0.0 && uz >= 0.0;
		inBox = inBox && ux <= 1.01 && uy <= 1.01 && uz <= 1.01;

		if (!inBox)
		{
			stringstream ss;
			ss << "encountered point outside bounding box." << endl;
			ss << "box.min: " << min.toString() << endl;
			ss << "box.max: " << max.toString() << endl;
			ss << "point: " << Vector3(x, y, z).toString() << endl;
			ss << "PotreeConverter requires a valid bounding box to operate." << endl;
			ss << "Please try to repair the bounding box, e.g. using lasinfo with the -repair_bb argument." << endl;
			logger::ERROR(ss.str());

			exit(123);
		}

		int64_t ix = int64_t(std::min(dGridSize * ux, dGridSize - 1.0));
		int64_t iy = int64_t(std::min(dGridSize * uy, dGridSize - 1.0));
		int64_t iz = int64_t(std::min(dGridSize * uz, dGridSize - 1.0));

		int64_t index = ix + iy * gridSize + iz * gridSize * gridSize;

		grid[index]++;
	};

	// grid contains index of node in nodes
	struct NodeLUT
	{
		int64_t gridSize;
		vector<int> grid;
	};

	CountResult countPointsInCells(vector<Source> sources, Vector3 min, Vector3 max, int64_t gridSize, State &state, Attributes &outputAttributes, Monitor *monitor, Options options)
	{

		cout << endl;
		cout << "=======================================" << endl;
		cout << "=== COUNTING                           " << endl;
		cout << "=======================================" << endl;

		auto tStart = now();

		// Vector3 size = max - min;

		vector<std::atomic_int32_t> grid(gridSize * gridSize * gridSize);
		vector<std::shared_ptr<ThreadSafeCluster>> gridCluster(gridSize * gridSize * gridSize);
		for (int i = 0; i < gridCluster.size(); ++i)
		{
			gridCluster[i] = std::make_shared<ThreadSafeCluster>();
		}

		struct Task
		{
			string path;
			int64_t totalPoints = 0;
			int64_t firstPoint;
			int64_t firstByte;
			int64_t numBytes;
			int64_t numPoints;
			int64_t bpp;
			Vector3 scale;
			Vector3 offset;
			Vector3 min;
			Vector3 max;
			Source source;
		};

		auto processor = [gridSize, &grid, tStart, &state, &outputAttributes, monitor, options, &gridCluster](shared_ptr<Task> task)
		{
			string path = task->path;
			int64_t start = task->firstByte;
			int64_t numBytes = task->numBytes;
			int64_t numToRead = task->numPoints;
			int64_t bpp = task->bpp;
			// Vector3 scale = task->scale;
			// Vector3 offset = task->offset;
			Vector3 min = task->min;
			Vector3 max = task->max;
			Source source = task->source;

			stringstream ss;
			ss << "counting " << fs::path(task->path).filename().string()
			   << ", first point: " << formatNumber(task->firstPoint)
			   << ", num points: " << formatNumber(task->numPoints);
			// cout << ss.str();
			// monitor->print("counter message", ss.str());

			logger::INFO(ss.str());
			// thread_local unique_ptr<void, void(*)(void*)> buffer(nullptr, free);
			// thread_local int64_t bufferSize = -1;

			// { // sanity checks
			// 	if(numBytes < 0){
			// 		logger::ERROR("invalid malloc size: " + formatNumber(numBytes));
			// 	}
			// }

			// if (bufferSize < numBytes){
			// 	buffer.reset(malloc(numBytes));
			// 	bufferSize = numBytes;
			// }

			double cubeSize = (max - min).max();
			Vector3 size = {cubeSize, cubeSize, cubeSize};
			max = min + cubeSize;

			double dGridSize = double(gridSize);

			double coordinates[3];

			auto posScale = outputAttributes.posScale;
			auto posOffset = outputAttributes.posOffset;

			string extension = std::filesystem::path(path).extension().string();
			if (extension == ".ply")
			{
				std::ifstream file(path, std::ios::binary);
				if (!file)
				{
					std::cerr << "Error: Unable to open file " << path << std::endl;
					exit(1234);
				}
				for (int i = 0; i < numToRead; i++)
				{
					int64_t currentPos = start + i * bpp;
					file.seekg(currentPos, std::ios::beg);
					char pointBuffer[12];
					file.read(pointBuffer, 12);
					if (!file)
						break;
					float x, y, z;
					memcpy(&x, &pointBuffer[0], sizeof(float));
					memcpy(&y, &pointBuffer[4], sizeof(float));
					memcpy(&z, &pointBuffer[8], sizeof(float));
					coordinates[0] = static_cast<double>(x);
					coordinates[1] = static_cast<double>(y);
					coordinates[2] = static_cast<double>(z);
					if (options.hasCluster && outputAttributes.get("clusterLabel"))
					{
						int label;
						char clusterLabel[4];
						file.seekg(currentPos + 248, std::ios::beg);
						file.read(clusterLabel, 4);
						float l;
						memcpy(&l, &clusterLabel[0], sizeof(float));
						label = static_cast<int>(l);
						countGridPointNumber(coordinates, grid, dGridSize, posOffset, posScale, min, max, size, gridCluster, label);
					}
					else
					{
						countGridPointNumber(coordinates, grid, dGridSize, posOffset, posScale, min, max, size);
					}
				}
				file.close();
			}
			else
			{
				laszip_POINTER laszip_reader;
				{
					laszip_BOOL is_compressed = iEndsWith(path, ".laz") ? 1 : 0;
					laszip_BOOL request_reader = 1;
					laszip_create(&laszip_reader);
					laszip_request_compatibility_mode(laszip_reader, request_reader);
					laszip_open_reader(laszip_reader, path.c_str(), &is_compressed);
					laszip_seek_point(laszip_reader, task->firstPoint);
				}

				for (int i = 0; i < numToRead; i++)
				{
					int64_t pointOffset = i * bpp;
					laszip_read_point(laszip_reader);
					laszip_get_coordinates(laszip_reader, coordinates);
					countGridPointNumber(coordinates, grid, dGridSize, posOffset, posScale, min, max, size);
				}
				laszip_close_reader(laszip_reader);
				laszip_destroy(laszip_reader);
			}
			static int64_t pointsProcessed = 0;
			pointsProcessed += task->numPoints;

			state.name = "COUNTING";
			state.pointsProcessed = pointsProcessed;
			state.duration = now() - tStart;

			// cout << ("end: " + formatNumber(dbgCurr)) << endl;
		};

		TaskPool<Task> pool(numChunkerThreads, processor);

		auto tStartTaskAssembly = now();

		for (auto source : sources)
		{
			// auto parallel = std::execution::par;
			// for_each(parallel, paths.begin(), paths.end(), [&mtx, &sources](string path) {
			int64_t bpp;
			int64_t numPoints;
			int64_t pointsLeft;
			int64_t batchSize = 1'000'000;
			int64_t numRead = 0;
			int64_t firstByte;
			int64_t numBytes;
			string extension = std::filesystem::path(source.path).extension().string();
			if (extension == ".ply")
			{
				PlyHeader header = std::get<PlyHeader>(source.header);
				bpp = header.byteLength;
				numPoints = header.numPoints;
				firstByte = header.headerEndIndex;
				pointsLeft = numPoints;
			}
			else
			{
				laszip_POINTER laszip_reader;
				laszip_header *header;
				// laszip_point* point;
				{
					laszip_create(&laszip_reader);

					laszip_BOOL request_reader = 1;
					laszip_BOOL is_compressed = iEndsWith(source.path, ".laz") ? 1 : 0;

					laszip_request_compatibility_mode(laszip_reader, request_reader);
					laszip_open_reader(laszip_reader, source.path.c_str(), &is_compressed);
					laszip_get_header_pointer(laszip_reader, &header);
				}
				bpp = header->point_data_record_length;
				numPoints = std::max(uint64_t(header->number_of_point_records), header->extended_number_of_point_records);
				pointsLeft = numPoints;
				firstByte = header->offset_to_point_data;
				laszip_close_reader(laszip_reader);
				laszip_destroy(laszip_reader);
			}

			while (pointsLeft > 0)
			{
				int64_t numToRead;
				if (pointsLeft < batchSize)
				{
					numToRead = pointsLeft;
					pointsLeft = 0;
				}
				else
				{
					numToRead = batchSize;
					pointsLeft = pointsLeft - batchSize;
				}
				int64_t numBytes = numToRead * bpp;
				auto task = make_shared<Task>();
				task->path = source.path;
				task->totalPoints = numPoints;
				task->firstPoint = numRead;
				task->firstByte = firstByte + numRead * bpp;
				task->numBytes = numBytes;
				task->numPoints = numToRead;
				task->bpp = bpp;
				// task->scale = { header->x_scale_factor, header->y_scale_factor, header->z_scale_factor };
				// task->offset = { header->x_offset, header->y_offset, header->z_offset };
				task->min = min;
				task->max = max;
				task->source = source;
				pool.addTask(task);

				numRead += batchSize;
			}
		}

		printElapsedTime("tStartTaskAssembly", tStartTaskAssembly);

		pool.waitTillEmpty();
		pool.close();

		printElapsedTime("countPointsInCells", tStart);

		double duration = now() - tStart;

		cout << "finished counting in " << formatNumber(duration) << "s" << endl;
		cout << "=======================================" << endl;

		{
			double duration = now() - tStart;
			state.values["duration(chunking-count)"] = formatNumber(duration, 3);
		}
		CountResult result;
		result.grid = std::move(grid);
		result.gridCluster = std::move(gridCluster);
		return result;
	}

	void addBuckets(string targetDir, vector<shared_ptr<Buffer>> &newBuckets)
	{

		for (int nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++)
		{

			if (newBuckets[nodeIndex]->size == 0)
			{
				continue;
			}

			auto &node = nodes[nodeIndex];
			string path = targetDir + "/chunks/" + node.id + ".bin";
			auto buffer = newBuckets[nodeIndex];

			writer->write(path, buffer);
		}
	}

	vector<function<void(int64_t)>> createAttributeHandlers(laszip_header *header, uint8_t *data, laszip_point *point, Attributes &inputAttributes, Attributes &outputAttributes)
	{

		vector<function<void(int64_t)>> handlers;
		int attributeOffset = 0;

		// reset min/max, we're writing the values to a per-thread copy anyway
		for (auto &attribute : outputAttributes.list)
		{
			attribute.min = {Infinity, Infinity, Infinity};
			attribute.max = {-Infinity, -Infinity, -Infinity};
		}

		{ // STANDARD LAS ATTRIBUTES

			int offsetRGB = outputAttributes.getOffset("rgb");
			Attribute *attributeRGB = outputAttributes.get("rgb");
			auto rgb = [data, point, header, offsetRGB, attributeRGB](int64_t offset)
			{
				if (offsetRGB >= 0)
				{

					uint16_t rgb[] = {0, 0, 0};
					memcpy(rgb, &point->rgb, 6);
					memcpy(data + offset + offsetRGB, rgb, 6);

					attributeRGB->min.x = std::min(attributeRGB->min.x, double(rgb[0]));
					attributeRGB->min.y = std::min(attributeRGB->min.y, double(rgb[1]));
					attributeRGB->min.z = std::min(attributeRGB->min.z, double(rgb[2]));

					attributeRGB->max.x = std::max(attributeRGB->max.x, double(rgb[0]));
					attributeRGB->max.y = std::max(attributeRGB->max.y, double(rgb[1]));
					attributeRGB->max.z = std::max(attributeRGB->max.z, double(rgb[2]));
				}
			};

			int offsetIntensity = outputAttributes.getOffset("intensity");
			Attribute *attributeIntensity = outputAttributes.get("intensity");
			auto intensity = [data, point, header, offsetIntensity, attributeIntensity](int64_t offset)
			{
				memcpy(data + offset + offsetIntensity, &point->intensity, 2);

				attributeIntensity->min.x = std::min(attributeIntensity->min.x, double(point->intensity));
				attributeIntensity->max.x = std::max(attributeIntensity->max.x, double(point->intensity));
			};

			int offsetReturnNumber = outputAttributes.getOffset("return number");
			Attribute *attributeReturnNumber = outputAttributes.get("return number");
			auto returnNumber = [data, point, header, offsetReturnNumber, attributeReturnNumber](int64_t offset)
			{
				uint8_t value = point->return_number;

				memcpy(data + offset + offsetReturnNumber, &value, 1);

				attributeReturnNumber->min.x = std::min(attributeReturnNumber->min.x, double(value));
				attributeReturnNumber->max.x = std::max(attributeReturnNumber->max.x, double(value));
			};

			int offsetNumberOfReturns = outputAttributes.getOffset("number of returns");
			Attribute *attributeNumberOfReturns = outputAttributes.get("number of returns");
			auto numberOfReturns = [data, point, header, offsetNumberOfReturns, attributeNumberOfReturns](int64_t offset)
			{
				uint8_t value = point->number_of_returns;

				memcpy(data + offset + offsetNumberOfReturns, &value, 1);

				attributeNumberOfReturns->min.x = std::min(attributeNumberOfReturns->min.x, double(value));
				attributeNumberOfReturns->max.x = std::max(attributeNumberOfReturns->max.x, double(value));
			};

			int offsetScanAngleRank = outputAttributes.getOffset("scan angle rank");
			Attribute *attributeScanAngleRank = outputAttributes.get("scan angle rank");
			auto scanAngleRank = [data, point, header, offsetScanAngleRank, attributeScanAngleRank](int64_t offset)
			{
				memcpy(data + offset + offsetScanAngleRank, &point->scan_angle_rank, 1);

				attributeScanAngleRank->min.x = std::min(attributeScanAngleRank->min.x, double(point->scan_angle_rank));
				attributeScanAngleRank->max.x = std::max(attributeScanAngleRank->max.x, double(point->scan_angle_rank));
			};

			int offsetScanAngle = outputAttributes.getOffset("scan angle");
			Attribute *attributeScanAngle = outputAttributes.get("scan angle");
			auto scanAngle = [data, point, header, offsetScanAngle, attributeScanAngle](int64_t offset)
			{
				memcpy(data + offset + offsetScanAngle, &point->extended_scan_angle, 2);

				attributeScanAngle->min.x = std::min(attributeScanAngle->min.x, double(point->extended_scan_angle));
				attributeScanAngle->max.x = std::max(attributeScanAngle->max.x, double(point->extended_scan_angle));
			};

			int offsetUserData = outputAttributes.getOffset("user data");
			Attribute *attributeUserData = outputAttributes.get("user data");
			auto userData = [data, point, header, offsetUserData, attributeUserData](int64_t offset)
			{
				memcpy(data + offset + offsetUserData, &point->user_data, 1);

				attributeUserData->min.x = std::min(attributeUserData->min.x, double(point->user_data));
				attributeUserData->max.x = std::max(attributeUserData->max.x, double(point->user_data));
			};

			int offsetClassification = outputAttributes.getOffset("classification");
			Attribute *attributeClassification = outputAttributes.get("classification");
			auto classification = [data, point, header, offsetClassification, attributeClassification](int64_t offset)
			{
				uint8_t value = 0;
				if (point->extended_classification > 31)
				{
					value = point->extended_classification;
				}
				else
				{
					value = point->classification;
				}

				data[offset + offsetClassification] = value;
				attributeClassification->histogram[value]++;

				attributeClassification->min.x = std::min(attributeClassification->min.x, double(value));
				attributeClassification->max.x = std::max(attributeClassification->max.x, double(value));
			};

			int offsetSourceId = outputAttributes.getOffset("point source id");
			Attribute *attributePointSourceId = outputAttributes.get("point source id");
			auto pointSourceId = [data, point, header, offsetSourceId, attributePointSourceId](int64_t offset)
			{
				memcpy(data + offset + offsetSourceId, &point->point_source_ID, 2);

				attributePointSourceId->min.x = std::min(attributePointSourceId->min.x, double(point->point_source_ID));
				attributePointSourceId->max.x = std::max(attributePointSourceId->max.x, double(point->point_source_ID));
			};

			int offsetGpsTime = outputAttributes.getOffset("gps-time");
			Attribute *attributeGpsTime = outputAttributes.get("gps-time");
			auto gpsTime = [data, point, header, offsetGpsTime, attributeGpsTime](int64_t offset)
			{
				memcpy(data + offset + offsetGpsTime, &point->gps_time, 8);

				attributeGpsTime->min.x = std::min(attributeGpsTime->min.x, point->gps_time);
				attributeGpsTime->max.x = std::max(attributeGpsTime->max.x, point->gps_time);
			};

			int offsetClassificationFlags = outputAttributes.getOffset("classification flags");
			Attribute *attributeClassificationFlags = outputAttributes.get("classification flags");
			auto classificationFlags = [data, point, header, offsetClassificationFlags, attributeClassificationFlags](int64_t offset)
			{
				uint8_t value = point->extended_classification_flags;

				memcpy(data + offset + offsetClassificationFlags, &value, 1);

				attributeClassificationFlags->min.x = std::min(attributeClassificationFlags->min.x, double(point->extended_classification_flags));
				attributeClassificationFlags->max.x = std::max(attributeClassificationFlags->max.x, double(point->extended_classification_flags));
			};

			unordered_map<string, function<void(int64_t)>> mapping = {
				{"rgb", rgb},
				{"intensity", intensity},
				{"return number", returnNumber},
				{"number of returns", numberOfReturns},
				{"classification", classification},
				{"scan angle rank", scanAngleRank},
				{"scan angle", scanAngle},
				{"user data", userData},
				{"point source id", pointSourceId},
				{"gps-time", gpsTime},
				{"classification flags", classificationFlags},
			};

			for (auto &attribute : inputAttributes.list)
			{

				attributeOffset += attribute.size;

				if (attribute.name == "position")
				{
					continue;
				}

				bool standardMappingExists = mapping.find(attribute.name) != mapping.end();
				bool isIncludedInOutput = outputAttributes.get(attribute.name) != nullptr;
				if (standardMappingExists && isIncludedInOutput)
				{
					handlers.push_back(mapping[attribute.name]);
				}
			}
		}

		{ // EXTRA ATTRIBUTES

			// mapping from las format to index of first extra attribute
			// +1 for all formats with returns, which is split into return number and number of returns
			unordered_map<int, int> formatToExtraIndex = {
				{0, 8},
				{1, 9},
				{2, 9},
				{3, 10},
				{4, 14},
				{5, 15},
				{6, 10},
				{7, 11},
			};

			bool noMapping = formatToExtraIndex.find(header->point_data_format) == formatToExtraIndex.end();
			if (noMapping)
			{
				string msg = "ERROR: las format not supported: " + formatNumber(header->point_data_format) + "\n";
				cout << msg;

				exit(123);
			}

			// handle extra bytes individually to compute per-attribute information
			int firstExtraIndex = formatToExtraIndex[header->point_data_format];
			int sourceOffset = 0;

			int attributeOffset = 0;
			for (int i = 0; i < firstExtraIndex; i++)
			{
				attributeOffset += inputAttributes.list[i].size;
			}

			for (int i = firstExtraIndex; i < inputAttributes.list.size(); i++)
			{
				Attribute &inputAttribute = inputAttributes.list[i];
				Attribute *attribute = outputAttributes.get(inputAttribute.name);
				int targetOffset = outputAttributes.getOffset(inputAttribute.name);

				int attributeSize = inputAttribute.size;

				if (attribute != nullptr)
				{
					auto handleAttribute = [data, point, header, attributeSize, attributeOffset, sourceOffset, attribute](int64_t offset)
					{
						memcpy(data + offset + attributeOffset, point->extra_bytes + sourceOffset, attributeSize);

						std::function<double(uint8_t *)> f;

						// TODO: shouldn't use DOUBLE as a unifying type
						// it won't work with uint64_t and int64_t
						if (attribute->type == AttributeType::INT8)
						{
							f = asDouble<int8_t>;
						}
						else if (attribute->type == AttributeType::INT16)
						{
							f = asDouble<int16_t>;
						}
						else if (attribute->type == AttributeType::INT32)
						{
							f = asDouble<int32_t>;
						}
						else if (attribute->type == AttributeType::INT64)
						{
							f = asDouble<int64_t>;
						}
						else if (attribute->type == AttributeType::UINT8)
						{
							f = asDouble<uint8_t>;
						}
						else if (attribute->type == AttributeType::UINT16)
						{
							f = asDouble<uint16_t>;
						}
						else if (attribute->type == AttributeType::UINT32)
						{
							f = asDouble<uint32_t>;
						}
						else if (attribute->type == AttributeType::UINT64)
						{
							f = asDouble<uint64_t>;
						}
						else if (attribute->type == AttributeType::FLOAT)
						{
							f = asDouble<float>;
						}
						else if (attribute->type == AttributeType::DOUBLE)
						{
							f = asDouble<double>;
						}

						if (attribute->numElements == 1)
						{
							double x = f(point->extra_bytes + sourceOffset);

							attribute->min.x = std::min(attribute->min.x, x);
							attribute->max.x = std::max(attribute->max.x, x);
						}
						else if (attribute->numElements == 2)
						{
							double x = f(point->extra_bytes + sourceOffset + 0 * attribute->elementSize);
							double y = f(point->extra_bytes + sourceOffset + 1 * attribute->elementSize);

							attribute->min.x = std::min(attribute->min.x, x);
							attribute->min.y = std::min(attribute->min.y, y);
							attribute->max.x = std::max(attribute->max.x, x);
							attribute->max.y = std::max(attribute->max.y, y);
						}
						else if (attribute->numElements == 3)
						{
							double x = f(point->extra_bytes + sourceOffset + 0 * attribute->elementSize);
							double y = f(point->extra_bytes + sourceOffset + 1 * attribute->elementSize);
							double z = f(point->extra_bytes + sourceOffset + 2 * attribute->elementSize);

							attribute->min.x = std::min(attribute->min.x, x);
							attribute->min.y = std::min(attribute->min.y, y);
							attribute->min.z = std::min(attribute->min.z, z);
							attribute->max.x = std::max(attribute->max.x, x);
							attribute->max.y = std::max(attribute->max.y, y);
							attribute->max.z = std::max(attribute->max.z, z);
						}
					};

					handlers.push_back(handleAttribute);
					attributeOffset += attribute->size;
				}

				sourceOffset += inputAttribute.size;
			}
		}

		return handlers;
	}

	void distributePoints(vector<Source> sources, Vector3 min, Vector3 max, string targetDir, NodeLUT &lut, State &state, Attributes &outputAttributes, Monitor *monitor, Options options)
	{

		cout << endl;
		cout << "=======================================" << endl;
		cout << "=== CREATING CHUNKS                    " << endl;
		cout << "=======================================" << endl;

		auto tStart = now();

		state.pointsProcessed = 0;
		state.bytesProcessed = 0;
		state.duration = 0;

		writer = new ConcurrentWriter(numFlushThreads, state);

		printElapsedTime("distributePoints0", tStart);

		vector<std::atomic_int32_t> counters(nodes.size());

		struct Task
		{
			string path;
			int64_t maxBatchSize;
			int64_t batchSize;
			int64_t firstPoint;
			NodeLUT *lut;
			Vector3 scale;
			Vector3 offset;
			Vector3 min;
			Vector3 max;
			Attributes inputAttributes;
			Source source;
		};

		mutex mtx_push_point;

		printElapsedTime("distributePoints1", tStart);

		auto processor = [&mtx_push_point, &counters, targetDir, &state, tStart, &outputAttributes, options](shared_ptr<Task> task)
		{
			auto path = task->path;
			auto batchSize = task->batchSize;
			auto *lut = task->lut;
			Source source = task->source;
			string extension = std::filesystem::path(source.path).extension().string();
			int bpp;
			if (extension == ".ply")
			{
				PlyHeader header = std::get<PlyHeader>(source.header);
				bpp = header.byteLength;
			}
			else
			{
				bpp = outputAttributes.bytes;
			}
			auto numBytes = bpp * batchSize;
			Vector3 scale = task->scale;
			Vector3 offset = task->offset;
			Vector3 min = task->min;
			Vector3 max = task->max;
			Attributes inputAttributes = task->inputAttributes;

			auto gridSize = lut->gridSize;
			auto &grid = lut->grid;

			std::unique_ptr<uint8_t[]> buffer;
			int64_t bufferSize = -1;

			if (numBytes < 0)
			{
				std::cerr << "Invalid malloc size: " << numBytes << std::endl;
				return;
			}

			if (bufferSize < numBytes)
			{
				buffer = std::make_unique<uint8_t[]>(numBytes);
				bufferSize = numBytes;
			}
			uint8_t *data = buffer.get();
			// memset necessary if attribute handlers don't set all values.
			// previous handlers from input with different point formats
			// may have set the values before.
			memset(data, 0, bufferSize);

			writer->waitUntilMemoryBelow(2'000);

			int pointFormat = -1;
			// per-thread copy of outputAttributes to compute min/max in a thread-safe way
			// will be merged to global outputAttributes instance at the end of this function
			Attributes outputAttributesCopy = outputAttributes;

			for (auto &attribute : outputAttributesCopy.list)
			{
				if (attribute.name == "classification")
				{
					for (int i = 0; i < attribute.histogram.size(); i++)
					{
						attribute.histogram[i] = 0;
					}
				}
			}
			if (extension == ".ply")
			{
				PlyHeader header = std::get<PlyHeader>(source.header);
				auto aPosition = outputAttributesCopy.get("position");
				std::ifstream file(path, std::ios::binary);
				if (!file)
				{
					std::cerr << "Error: Unable to open file " << path << std::endl;
					exit(1234);
				}
				for (int i = 0; i < batchSize; i++)
				{
					int64_t offset = i * header.byteLength;
					int bytesRemain = header.byteLength - 12;
					int64_t currentPos = header.headerEndIndex + task->firstPoint * header.byteLength + i * header.byteLength;
					file.seekg(currentPos, std::ios::beg);
					char pointBuffer[12];
					file.read(pointBuffer, 12);
					file.seekg(currentPos + 12, std::ios::beg);
					if (!file)
						break;
					float x, y, z;
					memcpy(&x, &pointBuffer[0], sizeof(float));
					memcpy(&y, &pointBuffer[4], sizeof(float));
					memcpy(&z, &pointBuffer[8], sizeof(float));
					x = (x - outputAttributes.posOffset.x) / scale.x;
					y = (y - outputAttributes.posOffset.y) / scale.y;
					z = (z - outputAttributes.posOffset.z) / scale.z;
					memcpy(data + offset + 0, &x, 4);
					memcpy(data + offset + 4, &y, 4);
					memcpy(data + offset + 8, &z, 4);
					if (options.SH == 3)
					{
						char *attribute = static_cast<char *>(malloc(bytesRemain));
						file.read(attribute, bytesRemain);
						memcpy(data + offset + 12, attribute, bytesRemain);
					}
					else if (options.SH == 2)
					{
						char *attribute_1 = static_cast<char *>(malloc(120));
						file.read(attribute_1, 120);
						bytesRemain = bytesRemain - 120 - 96;
						char *attribute_2 = static_cast<char *>(malloc(bytesRemain));
						int offset = outputAttributes.getOffset("opacity");
						file.seekg(currentPos + offset, std::ios::beg);
						file.read(attribute_2, bytesRemain);
						memcpy(data + offset + 12, attribute_1, 120);
						memcpy(data + offset + 132, &attribute_2, bytesRemain);
					}
					else if (options.SH == 1)
					{
						char *attribute_1 = static_cast<char *>(malloc(60));
						file.read(attribute_1, 60);
						bytesRemain = bytesRemain - 60 - 156;
						char *attribute_2 = static_cast<char *>(malloc(bytesRemain));
						int offset = outputAttributes.getOffset("opacity");
						file.seekg(currentPos + offset, std::ios::beg);
						file.read(attribute_2, bytesRemain);
						memcpy(data + offset + 12, attribute_1, 60);
						memcpy(data + offset + 132, &attribute_2, bytesRemain);
					}
					else
					{
						char *attribute_1 = static_cast<char *>(malloc(24));
						file.read(attribute_1, 24);
						bytesRemain = bytesRemain - 24 - 192;
						char *attribute_2 = static_cast<char *>(malloc(bytesRemain));
						int offset = outputAttributes.getOffset("opacity");
						file.seekg(currentPos + offset, std::ios::beg);
						file.read(attribute_2, bytesRemain);
						memcpy(data + offset + 12, attribute_1, 24);
						memcpy(data + offset + 132, &attribute_2, bytesRemain);
					}
				}
				file.close();
			}
			else
			{
				laszip_POINTER laszip_reader;
				laszip_header *header;
				laszip_point *point;

				laszip_BOOL request_reader = 1;
				laszip_BOOL is_compressed = iEndsWith(path, ".laz") ? 1 : 0;

				laszip_create(&laszip_reader);
				laszip_request_compatibility_mode(laszip_reader, request_reader);
				laszip_open_reader(laszip_reader, path.c_str(), &is_compressed);
				laszip_get_header_pointer(laszip_reader, &header);
				laszip_get_point_pointer(laszip_reader, &point);

				laszip_seek_point(laszip_reader, task->firstPoint);

				auto attributeHandlers = createAttributeHandlers(header, data, point, inputAttributes, outputAttributesCopy);

				double coordinates[3];
				auto aPosition = outputAttributesCopy.get("position");

				for (int64_t i = 0; i < batchSize; i++)
				{
					laszip_read_point(laszip_reader);
					laszip_get_coordinates(laszip_reader, coordinates);

					int64_t offset = i * outputAttributes.bytes;

					{ // copy position
						double x = coordinates[0];
						double y = coordinates[1];
						double z = coordinates[2];

						int32_t X = int32_t((x - outputAttributes.posOffset.x) / scale.x);
						int32_t Y = int32_t((y - outputAttributes.posOffset.y) / scale.y);
						int32_t Z = int32_t((z - outputAttributes.posOffset.z) / scale.z);

						memcpy(data + offset + 0, &X, 4);
						memcpy(data + offset + 4, &Y, 4);
						memcpy(data + offset + 8, &Z, 4);

						aPosition->min.x = std::min(aPosition->min.x, x);
						aPosition->min.y = std::min(aPosition->min.y, y);
						aPosition->min.z = std::min(aPosition->min.z, z);

						aPosition->max.x = std::max(aPosition->max.x, x);
						aPosition->max.y = std::max(aPosition->max.y, y);
						aPosition->max.z = std::max(aPosition->max.z, z);
					}

					// copy other attributes
					for (auto &handler : attributeHandlers)
					{
						handler(offset);
					}
				}

				pointFormat = header->point_data_format;

				laszip_close_reader(laszip_reader);
				laszip_destroy(laszip_reader);
			}

			double cubeSize = (max - min).max();
			Vector3 size = {cubeSize, cubeSize, cubeSize};
			max = min + cubeSize;

			double dGridSize = double(gridSize);

			auto toIndex = [data, &outputAttributes, scale, gridSize, dGridSize, size, min, path](int64_t pointOffset)
			{
				int64_t index;
				string extension = std::filesystem::path(path).extension().string();
				if (extension == ".ply")
				{
					float_t *xyz = reinterpret_cast<float_t *>(&data[0] + pointOffset);

					float_t X = xyz[0];
					float_t Y = xyz[1];
					float_t Z = xyz[2];

					double ux = (double(X) * scale.x + outputAttributes.posOffset.x - min.x) / size.x;
					double uy = (double(Y) * scale.y + outputAttributes.posOffset.y - min.y) / size.y;
					double uz = (double(Z) * scale.z + outputAttributes.posOffset.z - min.z) / size.z;

					int64_t ix = int64_t(std::min(dGridSize * ux, dGridSize - 1.0));
					int64_t iy = int64_t(std::min(dGridSize * uy, dGridSize - 1.0));
					int64_t iz = int64_t(std::min(dGridSize * uz, dGridSize - 1.0));

					index = ix + iy * gridSize + iz * gridSize * gridSize;
				}
				else
				{
					int32_t *xyz = reinterpret_cast<int32_t *>(&data[0] + pointOffset);

					int32_t X = xyz[0];
					int32_t Y = xyz[1];
					int32_t Z = xyz[2];

					double ux = (double(X) * scale.x + outputAttributes.posOffset.x - min.x) / size.x;
					double uy = (double(Y) * scale.y + outputAttributes.posOffset.y - min.y) / size.y;
					double uz = (double(Z) * scale.z + outputAttributes.posOffset.z - min.z) / size.z;

					int64_t ix = int64_t(std::min(dGridSize * ux, dGridSize - 1.0));
					int64_t iy = int64_t(std::min(dGridSize * uy, dGridSize - 1.0));
					int64_t iz = int64_t(std::min(dGridSize * uz, dGridSize - 1.0));

					index = ix + iy * gridSize + iz * gridSize * gridSize;
				}

				return index;
			};

			// COUNT POINTS PER BUCKET
			vector<int64_t> counts(nodes.size(), 0);
			for (int64_t i = 0; i < batchSize; i++)
			{
				auto index = toIndex(i * bpp);

				auto nodeIndex = grid[index];

				// ERROR
				if (nodeIndex == -1)
				{

					int32_t *xyz = reinterpret_cast<int32_t *>(&data[0] + i * bpp);

					auto x = xyz[0];
					auto y = xyz[1];
					auto z = xyz[2];

					double ux = (xyz[0] * scale.x + outputAttributes.posOffset.x) / size.x;
					double uy = (xyz[1] * scale.y + outputAttributes.posOffset.y) / size.y;
					double uz = (xyz[2] * scale.z + outputAttributes.posOffset.z) / size.z;

					int64_t ix = int64_t(std::min(dGridSize * ux, dGridSize - 1.0));
					int64_t iy = int64_t(std::min(dGridSize * uy, dGridSize - 1.0));
					int64_t iz = int64_t(std::min(dGridSize * uz, dGridSize - 1.0));

					int64_t index = ix + iy * gridSize + iz * gridSize * gridSize;

					stringstream ss;
					ss << "point to node lookup failed, no node found." << endl;
					ss << "point: " << formatNumber(x, 3) << ", " << formatNumber(y, 3) << ", " << formatNumber(z, 3) << endl;
					ss << "3d grid index: " << ix << ", " << iy << ", " << iz << endl;
					ss << "1d grid index: " << index << endl;

					logger::ERROR(ss.str());
					exit(123);
				}

				counts[nodeIndex]++;
			}

			// ALLOCATE BUCKETS
			vector<shared_ptr<Buffer>> buckets(nodes.size(), nullptr);
			for (int i = 0; i < nodes.size(); i++)
			{
				int64_t numPoints = counts[i];
				int64_t bytes = numPoints * outputAttributes.bytes;
				buckets[i] = make_shared<Buffer>(bytes);
			}

			// ADD POINTS TO BUCKETS
			shared_ptr<Buffer> previousBucket = nullptr;
			int64_t previousNodeIndex = -1;
			for (int64_t i = 0; i < batchSize; i++)
			{
				int64_t pointOffset = i * bpp;

				auto index = toIndex(pointOffset);

				auto nodeIndex = grid[index];
				auto &node = nodes[nodeIndex];

				if (nodeIndex == previousNodeIndex)
				{
					previousBucket->write(&data[0] + pointOffset, outputAttributes.bytes);
				}
				else
				{
					previousBucket = buckets[nodeIndex];
					previousNodeIndex = nodeIndex;
					previousBucket->write(&data[0] + pointOffset, outputAttributes.bytes);
				}
			}

			state.pointsProcessed += batchSize;
			state.bytesProcessed += numBytes;
			state.duration = now() - tStart;

			auto tAddBuckets = now();
			addBuckets(targetDir, buckets);

			// merge attribute metadata of this batch into global attribute metadata
			for (int i = 0; i < outputAttributesCopy.list.size(); i++)
			{
				Attribute &source = outputAttributesCopy.list[i];
				Attribute &target = outputAttributes.list[i];

				lock_guard<mutex> lock(mtx_attributes);
				target.min.x = std::min(target.min.x, source.min.x);
				target.min.y = std::min(target.min.y, source.min.y);
				target.min.z = std::min(target.min.z, source.min.z);

				target.max.x = std::max(target.max.x, source.max.x);
				target.max.y = std::max(target.max.y, source.max.y);
				target.max.z = std::max(target.max.z, source.max.z);

				// target.mask = target.mask | source.mask;

				for (int j = 0; j < target.histogram.size(); j++)
				{
					target.histogram[j] = target.histogram[j] + source.histogram[j];
				}
			}
		};

		TaskPool<Task> pool(numChunkerThreads, processor);

		for (auto source : sources)
		{
			string extension = std::filesystem::path(source.path).extension().string();
			int64_t numPoints;
			int64_t pointsLeft;
			int64_t maxBatchSize = 1'000'000;
			int64_t numRead = 0;
			vector<Source> tmpSources = {source};
			Attributes inputAttributes = computeOutputAttributes(tmpSources, {}, options);
			if (extension == ".ply")
			{
				auto header = std::get<PlyHeader>(source.header);
				numPoints = header.numPoints;
				pointsLeft = numPoints;
			}
			else
			{
				laszip_POINTER laszip_reader;
				laszip_header *header;
				{
					laszip_BOOL request_reader = 1;
					laszip_BOOL is_compressed = iEndsWith(source.path, ".laz") ? 1 : 0;

					laszip_create(&laszip_reader);
					laszip_request_compatibility_mode(laszip_reader, request_reader);
					laszip_open_reader(laszip_reader, source.path.c_str(), &is_compressed);
					laszip_get_header_pointer(laszip_reader, &header);
				}
				numPoints = std::max(uint64_t(header->number_of_point_records), header->extended_number_of_point_records);
				pointsLeft = numPoints;
				laszip_close_reader(laszip_reader);
				laszip_destroy(laszip_reader);
			}
			while (pointsLeft > 0)
			{
				int64_t numToRead;
				if (pointsLeft < maxBatchSize)
				{
					numToRead = pointsLeft;
					pointsLeft = 0;
				}
				else
				{
					numToRead = maxBatchSize;
					pointsLeft = pointsLeft - maxBatchSize;
				}
				auto task = make_shared<Task>();
				task->maxBatchSize = maxBatchSize;
				task->batchSize = numToRead;
				task->lut = &lut;
				task->firstPoint = numRead;
				task->path = source.path;
				// task->scale = { header->x_scale_factor, header->y_scale_factor, header->z_scale_factor };
				task->scale = outputAttributes.posScale;
				task->offset = outputAttributes.posOffset;
				task->min = min;
				task->max = max;
				task->inputAttributes = inputAttributes;
				task->source = source;

				pool.addTask(task);

				numRead += numToRead;
			}
		}

		pool.close();
		writer->join();

		delete writer;

		auto duration = now() - tStart;
		cout << "finished creating chunks in " << formatNumber(duration) << "s" << endl;
		cout << "=======================================" << endl;
	}

	void writeMetadata(string path, Vector3 min, Vector3 max, Attributes &attributes)
	{
		json js;

		js["min"] = {min.x, min.y, min.z};
		js["max"] = {max.x, max.y, max.z};

		js["attributes"] = {};
		for (auto attribute : attributes.list)
		{

			json jsAttribute;
			jsAttribute["name"] = attribute.name;
			jsAttribute["size"] = attribute.size;
			jsAttribute["numElements"] = attribute.numElements;
			jsAttribute["elementSize"] = attribute.elementSize;
			jsAttribute["description"] = attribute.description;
			jsAttribute["type"] = getAttributeTypename(attribute.type);

			if (attribute.numElements == 1)
			{
				jsAttribute["min"] = vector<double>{attribute.min.x};
				jsAttribute["max"] = vector<double>{attribute.max.x};
				jsAttribute["scale"] = vector<double>{attribute.scale.x};
				jsAttribute["offset"] = vector<double>{attribute.offset.x};
			}
			else if (attribute.numElements == 2)
			{
				jsAttribute["min"] = vector<double>{attribute.min.x, attribute.min.y};
				jsAttribute["max"] = vector<double>{attribute.max.x, attribute.max.y};
				jsAttribute["scale"] = vector<double>{attribute.scale.x, attribute.scale.y};
				jsAttribute["offset"] = vector<double>{attribute.offset.x, attribute.offset.y};
			}
			else if (attribute.numElements == 3)
			{
				jsAttribute["min"] = vector<double>{attribute.min.x, attribute.min.y, attribute.min.z};
				jsAttribute["max"] = vector<double>{attribute.max.x, attribute.max.y, attribute.max.z};
				jsAttribute["scale"] = vector<double>{attribute.scale.x, attribute.scale.y, attribute.scale.z};
				jsAttribute["offset"] = vector<double>{attribute.offset.x, attribute.offset.y, attribute.offset.z};
			}

			bool emptyHistogram = true;
			for (int i = 0; i < attribute.histogram.size(); i++)
			{
				if (attribute.histogram[i] != 0)
				{
					emptyHistogram = false;
				}
			}

			if (attribute.size == 1 && !emptyHistogram)
			{
				json jsHistogram = attribute.histogram;

				jsAttribute["histogram"] = jsHistogram;
			}

			js["attributes"].push_back(jsAttribute);
		}

		js["scale"] = vector<double>({attributes.posScale.x,
									  attributes.posScale.y,
									  attributes.posScale.z});

		js["offset"] = vector<double>({attributes.posOffset.x,
									   attributes.posOffset.y,
									   attributes.posOffset.z});
		// js["nodes"] = {};
		// for (auto node : nodes)
		// {
		// 	json nodeInfo;
		// 	nodeInfo["id"] = node.id;
		// 	nodeInfo["level"] = node.level;
		// 	nodeInfo["x"] = node.x;
		// 	nodeInfo["y"] = node.y;
		// 	nodeInfo["z"] = node.z;
		// 	nodeInfo["size"] = node.size;
		// 	nodeInfo["numPoints"] = node.numPoints;
		// 	nodeInfo["clusterPoints"] = {};
		// 	for (auto &pair : node.clusterMap)
		// 	{
		// 		json clusterInfo;
		// 		int key = pair.first;
		// 		int value = pair.second;
		// 		clusterInfo[std::to_string(key)] = value;
		// 		nodeInfo["clusterPoints"].push_back(clusterInfo);
		// 	}
		// 	js["nodes"].push_back(nodeInfo);
		// }

		string content = js.dump(4);

		writeFile(path, content);
	}

	void writeSampleData(string path, vector<std::unordered_map<int, double>> sampleRate)
	{
		json js;

		for (int i = 1; i < sampleRate.size() + 1; i++)
		{
			string level = std::to_string(i);
			js[level] = {};
			for (auto categoryItem : sampleRate[i - 1])
			{
				json categorySampleRate;
				int key = categoryItem.first;
				double value = categoryItem.second;
				categorySampleRate[std::to_string(key)] = value;
				js[level].push_back(categorySampleRate);
			}
		}
		string content = js.dump(4);

		writeFile(path, content);
	}

	//
	// XXX_high: variables of the higher/more detailed level of the pyramid that we're evaluating right now
	// XXX_low: one level lower than _high; the target of the "downsampling" operation
	//
	NodeLUT createLUT(vector<atomic_int32_t> &grid, int64_t gridSize)
	{
		auto tStart = now();

		auto for_xyz = [](int64_t gridSize, function<void(int64_t, int64_t, int64_t)> callback)
		{
			for (int x = 0; x < gridSize; x++)
			{
				for (int y = 0; y < gridSize; y++)
				{
					for (int z = 0; z < gridSize; z++)
					{
						callback(x, y, z);
					}
				}
			}
		};

		// atomic vectors are cumbersome, convert the highest level into a regular integer vector first.
		vector<int64_t> grid_high;
		grid_high.reserve(grid.size());
		for (auto &value : grid)
		{
			grid_high.push_back(value);
		}

		int64_t level_max = int64_t(log2(gridSize));

		// - evaluate counting grid in "image pyramid" fashion
		// - merge smaller cells into larger ones
		// - unmergeable cells are resulting chunks; push them to "nodes" array.
		for (int64_t level_low = level_max - 1; level_low >= 0; level_low--)
		{

			int64_t level_high = level_low + 1;

			int64_t gridSize_high = pow(2, level_high);
			int64_t gridSize_low = pow(2, level_low);

			vector<int64_t> grid_low(gridSize_low * gridSize_low * gridSize_low, 0);
			// grid_high

			// loop through all cells of the lower detail target grid, and for each cell through the 8 enclosed cells of the higher level grid
			for_xyz(gridSize_low, [&grid_low, &grid_high, gridSize_low, gridSize_high, level_low, level_high, level_max](int64_t x, int64_t y, int64_t z)
					{
						int64_t index_low = x + y * gridSize_low + z * gridSize_low * gridSize_low;

						int sum = 0;
						int max = 0;
						bool unmergeable = false;

						// loop through the 8 enclosed cells of the higher detailed grid
						for (int64_t j = 0; j < 8; j++)
						{
							int64_t ox = (j & 0b100) >> 2;
							int64_t oy = (j & 0b010) >> 1;
							int64_t oz = (j & 0b001) >> 0;

							int64_t nx = 2 * x + ox;
							int64_t ny = 2 * y + oy;
							int64_t nz = 2 * z + oz;

							int64_t index_high = nx + ny * gridSize_high + nz * gridSize_high * gridSize_high;

							auto value = grid_high[index_high];
							int intValue = static_cast<int>(value);

							if (value == -1)
							{
								unmergeable = true;
							}
							else
							{
								sum += intValue;
							}

							max = std::max(max, intValue);
						}

						if (unmergeable || (sum - 100000) > 0)
						{

							// finished chunks
							for (int64_t j = 0; j < 8; j++)
							{
								int64_t ox = (j & 0b100) >> 2;
								int64_t oy = (j & 0b010) >> 1;
								int64_t oz = (j & 0b001) >> 0;

								int64_t nx = 2 * x + ox;
								int64_t ny = 2 * y + oy;
								int64_t nz = 2 * z + oz;

								int64_t index_high = nx + ny * gridSize_high + nz * gridSize_high * gridSize_high;

								auto value = grid_high[index_high];

								if (value > 0)
								{
									string nodeID = toNodeID(level_high, gridSize_high, nx, ny, nz);
									Node node(nodeID, value);
									node.x = nx;
									node.y = ny;
									node.z = nz;
									node.size = pow(2, (level_max - level_high));

									nodes.push_back(node);
								}
							}

							// invalidate the field to show the parent that nothing can be merged with it
							grid_low[index_low] = -1;
						}
						else
						{
							grid_low[index_low] = sum;
						} });

			grid_high = grid_low;
		}

		// - create lookup table
		// - loop through nodes, add pointers to node/chunk for all enclosed cells in LUT.
		vector<int32_t> lut(gridSize * gridSize * gridSize, -1);
		for (int i = 0; i < nodes.size(); i++)
		{
			auto node = nodes[i];

			for_xyz(node.size, [node, &lut, gridSize, i](int64_t ox, int64_t oy, int64_t oz)
					{
				int64_t x = node.size * node.x + ox;
				int64_t y = node.size * node.y + oy;
				int64_t z = node.size * node.z + oz;
				int64_t index = x + y * gridSize + z * gridSize * gridSize;

				lut[index] = i; });
		}

		printElapsedTime("createLUT", tStart);

		return {gridSize, lut};
	}

	void calculateSampleRate(vector<std::unordered_map<int, double>> &sampleRate, int level_max,
							 vector<std::shared_ptr<ThreadSafeCluster>> gridCluster, std::unordered_map<int, chunker_countsort_laszip::ClusterInfo> clusterStatistics)
	{
		auto for_xyz = [](int64_t gridSize, function<void(int64_t, int64_t, int64_t)> callback)
		{
			for (int x = 0; x < gridSize; x++)
			{
				for (int y = 0; y < gridSize; y++)
				{
					for (int z = 0; z < gridSize; z++)
					{
						callback(x, y, z);
					}
				}
			}
		};
		vector<std::unordered_map<int, int>> grid_cluster_high;
		grid_cluster_high.reserve(gridCluster.size());
		for (auto &value : gridCluster)
		{
			grid_cluster_high.push_back(std::move(value->clusterMap));
		}
		for (int level_low = level_max - 1; level_low >= 0; level_low--)
		{
			int level_high = level_low + 1;
			int gridSizeHigh = pow(2, level_high);
			int gridSizeLow = pow(2, level_low);
			vector<std::unordered_map<int, int>> grid_cluster_low(gridSizeLow * gridSizeLow * gridSizeLow);

			// loop through all cells of the lower detail target grid, and for each cell through the 8 enclosed cells of the higher level grid
			for_xyz(gridSizeLow, [gridSizeLow, gridSizeHigh, level_low, level_high, level_max, &grid_cluster_high, &grid_cluster_low, &clusterStatistics, &sampleRate](int64_t x, int64_t y, int64_t z)
					{
						int64_t index_low = x + y * gridSizeLow + z * gridSizeLow * gridSizeLow;
						std::unordered_map<int, int> clusterSum;

						// loop through the 8 enclosed cells of the higher detailed grid
						for (int64_t j = 0; j < 8; j++)
						{
							float weightingSum = 0.0;
							int sum = 0;
							std::unordered_map<int, float> clusterWeighting;
							int64_t ox = (j & 0b100) >> 2;
							int64_t oy = (j & 0b010) >> 1;
							int64_t oz = (j & 0b001) >> 0;

							int64_t nx = 2 * x + ox;
							int64_t ny = 2 * y + oy;
							int64_t nz = 2 * z + oz;

							int64_t index_high = nx + ny * gridSizeHigh + nz * gridSizeHigh * gridSizeHigh;

							auto clusterValue = grid_cluster_high[index_high];

							for (auto& pair : clusterValue) {
								int key = pair.first;
								int value = pair.second;     // 值
								double dis = clusterStatistics[key].dis;
								double weighting = (1.0/dis)*value;
								clusterWeighting[key] = weighting;
								sum += value;
								weightingSum += weighting;
								clusterSum[key] += value;
							}
							for (auto& pair : clusterValue) {
								int key = pair.first;
								int value = pair.second;     // 值
								clusterWeighting[key] = (clusterWeighting[key]*sum)/(clusterValue[key]*weightingSum);
								if(clusterWeighting[key]>sampleRate[level_high - 1][key]){
									sampleRate[level_high - 1][key] = clusterWeighting[key];
								} 
							}
						}
							grid_cluster_low[index_low] = clusterSum; });
			grid_cluster_high = grid_cluster_low;
		}
	}

	void doChunking(vector<Source> sources, string targetDir, Vector3 min, Vector3 max, State &state,
					Attributes outputAttributes, Monitor *monitor, Options options, std::unordered_map<int, chunker_countsort_laszip::ClusterInfo> clusterInfoMap)
	{

		auto tStart = now();

		int64_t tmp = state.pointsTotal / 20;
		maxPointsPerChunk = std::min(tmp, int64_t(10'000'000));
		// cout << "maxPointsPerChunk: " << maxPointsPerChunk << endl;

		if (state.pointsTotal < 1'000'000)
		{
			gridSize = 128;
		}
		else if (state.pointsTotal < 5'000'000)
		{
			gridSize = 256;
		}
		else
		{
			gridSize = 512;
		}

		state.currentPass = 1;

		{ // prepare/clean target directories
			string dir = targetDir + "/chunks";
			fs::create_directories(dir);

			for (const auto &entry : std::filesystem::directory_iterator(dir))
			{
				std::filesystem::remove(entry);
			}
		}

		// COUNT
		auto countResult = countPointsInCells(sources, min, max, gridSize, state, outputAttributes, monitor, options);

		int level_max = log2(gridSize);
		vector<std::unordered_map<int, double>> sampleRate(level_max);

		calculateSampleRate(sampleRate, level_max, countResult.gridCluster,clusterInfoMap);

		{ // DISTIRBUTE
			auto tStartDistribute = now();
			auto lut = createLUT(countResult.grid, gridSize);

			state.currentPass = 2;
			distributePoints(sources, min, max, targetDir, lut, state, outputAttributes, monitor, options);

			{
				double duration = now() - tStartDistribute;
				state.values["duration(chunking-distribute)"] = formatNumber(duration, 3);
			}
		}

		string metadataPath = targetDir + "/chunks/metadata.json";
		double cubeSize = (max - min).max();
		Vector3 size = {cubeSize, cubeSize, cubeSize};
		max = min + cubeSize;

		writeMetadata(metadataPath, min, max, outputAttributes);

		string sampleDataPath =  targetDir + "/sampleData.json";
		writeSampleData(sampleDataPath, sampleRate);

		double duration = now() - tStart;
		state.values["duration(chunking-total)"] = formatNumber(duration, 3);
	}

}
