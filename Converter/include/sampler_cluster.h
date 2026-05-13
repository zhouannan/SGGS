
#pragma once

#include <execution>
#include <utility> // 包含 std::pair
#include "structures.h"
#include "Attributes.h"
#include "converter_utils.h"

struct SamplerCluster : public Sampler
{

	float sigmoid(float x)
	{
		return 1.0 / (1.0 + exp(-x)); // exp函数计算e的x次方
	};
	// subsample a local octree from bottom up
	void sample(Node *node, Attributes attributes, double baseSpacing,
				function<void(Node *)> onNodeCompleted,
				function<void(Node *)> onNodeDiscarded,
				Options options)
	{
		auto sampleData = options.sampleRate;
		auto levelSampleRate = options.levelSampleRate;
		function<void(Node *, function<void(Node *)>)> traversePost = [&traversePost](Node *node, function<void(Node *)> callback)
		{
			for (auto child : node->children)
			{

				if (child != nullptr && !child->sampled)
				{
					traversePost(child.get(), callback);
				}
			}

			callback(node);
		};

		int bytesPerPoint = attributes.bytes;

		traversePost(node, [bytesPerPoint, &onNodeCompleted, &onNodeDiscarded, attributes, options, &sampleData, &levelSampleRate,this](Node *node)
					 {
			node->sampled = true;

			int64_t numPoints = node->numPoints;

			int64_t gridSize = 128;

			bool isLeaf = node->isLeaf();
			if (isLeaf) {
				// shuffle?

				//
				// a not particularly efficient approach to shuffling:
				// 

				vector<int> indices(node->numPoints);
				for (int i = 0; i < node->numPoints; i++) {
					indices[i] = i;
				}

				unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

				shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));

				auto buffer = make_shared<Buffer>(node->points->size);

				for (int i = 0; i < node->numPoints; i++) {

					int64_t sourceOffset = i * attributes.bytes;
					int64_t targetOffset = indices[i] * attributes.bytes;

					memcpy(buffer->data_u8 + targetOffset, node->points->data_u8 + sourceOffset, attributes.bytes);

				}

				node->points = buffer;


				return false;
			}

			// =================================================================
			// SAMPLING
			// =================================================================
			//
			// first, check for each point whether it's accepted or rejected
			// save result in an array with one element for each point

			vector<vector<int8_t>> acceptedChildPointFlags;
			vector<int64_t> numRejectedPerChild;
			int64_t numAccepted = 0;
			for (int childIndex = 0; childIndex < 8; childIndex++) {
				auto child = node->children[childIndex];
				if (child == nullptr) {
					acceptedChildPointFlags.push_back({});
					numRejectedPerChild.push_back({});

					continue;
				}
				int level = child->level();
				if(level>8)level = 8;
				auto clusterSampleData = sampleData[level];
				auto sampleRate = levelSampleRate[level];
				vector<int8_t> acceptedFlags(child->numPoints, 0);
				unordered_map<int,std::vector<std::pair<int,float>>> weightMap;
				unordered_map<int,bool>topIndices;
				int64_t numRejected = 0;

				for (int i = 0; i < child->numPoints; i++) {
					int64_t pointOffset = i * attributes.bytes;
					float_t* opacityAndScale = reinterpret_cast<float_t*>(child->points->data_u8 + pointOffset + 216);
					float_t* Category = reinterpret_cast<float_t*>(child->points->data_u8 + pointOffset + 248);
					float opacity = sigmoid(opacityAndScale[0]);
					float scaleX = exp(opacityAndScale[1]);
					float scaleY = exp(opacityAndScale[2]);
					float scaleZ = exp(opacityAndScale[3]);
					int category = (int)Category[0];
					float maxInfluence = std::max({scaleX * scaleY, scaleX * scaleZ, scaleY * scaleZ});
					weightMap[category].push_back({i,maxInfluence});
				}
				for(auto weight:weightMap){
					int category = weight.first;
					double rate = clusterSampleData[category]*sampleRate;
					// double rate = 0.6;
					if(rate>1.0)rate=1.0;
					auto weightArray = weight.second;
					// 按照值进行排序
					std::sort(weightArray.begin(), weightArray.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
						return a.second > b.second; // 降序排序
					});
					// 创建随机数引擎
					// std::random_device rd;  // 用于生成种子
					// std::mt19937 g(rd());   // 生成随机数引擎
					// std::shuffle(weightArray.begin(), weightArray.end(), g);
					// 计算前 70% 的元素数量
					int numToKeep = static_cast<int>(weightArray.size() * rate);

					// 提取前 70% 的索引
					
					for (int i = 0; i < numToKeep; i++) {
						topIndices[weightArray[i].first] = true; // 添加索引
					}

				}
				for (int i = 0; i < child->numPoints; i++) {
					bool isAccepted;
					if (topIndices[i]) {
						isAccepted = true;
					} else {
						isAccepted = false;
					}

					if (isAccepted) {
						numAccepted++;
					} else {
						numRejected++;
					}

					acceptedFlags[i] = isAccepted ? 1 : 0;
				}
				acceptedChildPointFlags.push_back(acceptedFlags);
				numRejectedPerChild.push_back(numRejected);
			}

			auto accepted = make_shared<Buffer>(numAccepted * attributes.bytes);
			for (int childIndex = 0; childIndex < 8; childIndex++) {
				auto child = node->children[childIndex];

				if (child == nullptr) continue;

				auto numRejected = numRejectedPerChild[childIndex];
				auto& acceptedFlags = acceptedChildPointFlags[childIndex];
				auto rejected = make_shared<Buffer>(numRejected * attributes.bytes);

				for (int i = 0; i < child->numPoints; i++) {
					auto isAccepted = acceptedFlags[i];
					int64_t pointOffset = i * attributes.bytes;

					if (isAccepted) {
						accepted->write(child->points->data_u8 + pointOffset, attributes.bytes);
					} else {
						rejected->write(child->points->data_u8 + pointOffset, attributes.bytes);
					}
				}

				if (numRejected == 0 && child->isLeaf()) {
					onNodeDiscarded(child.get());

					node->children[childIndex] = nullptr;
				} if (numRejected > 0) {
					child->points = rejected;
					child->numPoints = numRejected;

					onNodeCompleted(child.get());
				} else if(numRejected == 0) {
					// the parent has taken all points from this child, 
					// so make this child an empty inner node.
					// Otherwise, the hierarchy file will claim that 
					// this node has points but because it doesn't have any,
					// decompressing the nonexistent point buffer fails
					// https://github.com/potree/potree/issues/1125
					child->points = nullptr;
					child->numPoints = 0;
					onNodeCompleted(child.get());
				}
			}

			node->points = accepted;
			node->numPoints = numAccepted;

			return true; 
			
			});
	}
};
