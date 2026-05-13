
#pragma once

#include "unsuck/unsuck.hpp"
#include "Vector3.h"
#include "Attributes.h"
#include <unordered_map>
#include <vector>


struct PlyTypeInfo {
	AttributeType type = AttributeType::UNDEFINED;
	int numElements = 0;
};

struct PlyHeader {
	unordered_map<string, string> types;
	unordered_map<string, int> bytes;
	int headerEndIndex;
	int byteLength;
	int64_t numPoints = 0;
	Vector3 min;
	Vector3 max;
};

struct PlyBox {
	Vector3 min={Infinity,Infinity,Infinity};
	Vector3 max={-Infinity,-Infinity,-Infinity};
};


PlyTypeInfo plyTypeInfo(int typeID);

int getSizeInBytes(const string& type);

PlyHeader loadPlyHeader(string path);

PlyBox getBoundingBox(string path,PlyHeader header);