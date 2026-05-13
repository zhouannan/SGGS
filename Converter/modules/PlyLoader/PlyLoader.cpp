
#include "PlyLoader/PlyLoader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <regex>

using namespace std;

unordered_map<string, string> TYPE_MAP = {
    {"double", "getFloat64"},
    {"int", "getInt32"},
    {"uint", "getUint32"},
    {"float", "getFloat32"},
    {"short", "getInt16"},
    {"ushort", "getUint16"},
    {"uchar", "getUint8"},
    {"char", "getInt8"}
};

unordered_map<string, int> SIZE_MAP = {
    {"double", 8},
    {"int", 4},
    {"uint", 4},
    {"float", 4},
    {"short", 2},
    {"ushort", 2},
    {"uchar", 1},
    {"char", 1}
};

int getSizeInBytes(const string& type) {
    if (type == "getFloat64") return 8;
    if (type == "getFloat32") return 4;
    if (type == "getInt32" || type == "getUint32") return 4;
    if (type == "getInt16" || type == "getUint16") return 2;
    if (type == "getUint8" || type == "getInt8") return 1;
    return 1; // Default to 1 byte for unknown types
}

PlyTypeInfo plyTypeInfo(int typeID) {

	unordered_map<int, AttributeType> mapping = {
		{0, AttributeType::UNDEFINED},
		{1, AttributeType::UINT8},
		{2, AttributeType::INT8},
		{3, AttributeType::UINT16},
		{4, AttributeType::INT16},
		{5, AttributeType::UINT32},
		{6, AttributeType::INT32},
		{7, AttributeType::UINT64},
		{8, AttributeType::INT64},
		{9, AttributeType::FLOAT},
		{10, AttributeType::DOUBLE},

		{11, AttributeType::UINT8},
		{12, AttributeType::INT8},
		{13, AttributeType::UINT16},
		{14, AttributeType::INT16},
		{15, AttributeType::UINT32},
		{16, AttributeType::INT32},
		{17, AttributeType::UINT64},
		{18, AttributeType::INT64},
		{19, AttributeType::FLOAT},
		{20, AttributeType::DOUBLE},

		{21, AttributeType::UINT8},
		{22, AttributeType::INT8},
		{23, AttributeType::UINT16},
		{24, AttributeType::INT16},
		{25, AttributeType::UINT32},
		{26, AttributeType::INT32},
		{27, AttributeType::UINT64},
		{28, AttributeType::INT64},
		{29, AttributeType::FLOAT},
		{30, AttributeType::DOUBLE},
	};

	if (mapping.find(typeID) != mapping.end()) {

		AttributeType type = mapping[typeID];

		int numElements = 0;
		if (typeID <= 10) {
			numElements = 1;
		} else if (typeID <= 20) {
			numElements = 2;
		} else if (typeID <= 30) {
			numElements = 3;
		}

		PlyTypeInfo info;
		info.type = type;
		info.numElements = numElements;

		return info;
	} else {
		cout << "ERROR: unkown extra attribute type: " << typeID << endl;
		exit(123);
	}


}


PlyHeader loadPlyHeader(string path) {
	PlyHeader result;
	ifstream file(path, ios::binary);

    if (!file.is_open()) {
        cerr << "Unable to open file!" << endl;
        return result;
    }

    // Read the first 10KB of the file for the header
    const size_t headerSize = 1024 * 10;
    vector<char> buffer(headerSize);
    file.read(buffer.data(), headerSize);
    string header(buffer.begin(), buffer.end());

    string headerEnd = "end_header\n";
    size_t headerEndIndex = header.find(headerEnd)+headerEnd.length();
    if (headerEndIndex == string::npos) {
        cerr << "Unable to read .ply file header" << endl;
        return result;
    }

    // Extract vertex count
	regex vertexRegex(R"(element\s+vertex\s+(\d+))");
    smatch match;
	int vertexCount;
    if (regex_search(header, match, vertexRegex)) {
        vertexCount = stoi(match[1].str());
        cout << "Vertex Count: " << vertexCount << endl;
    } else {
        cerr << "Vertex count not found in the header" << endl;
        return result;
    }
    // Parse properties and offsets
    unordered_map<string, string> types;
    unordered_map<string, int> bytes;
    int rowOffset = 0;

    stringstream headerStream(header.substr(0, headerEndIndex));
    string line;

	while (getline(headerStream, line)) {
        if (line.find("property ") == 0) {
            istringstream lineStream(line);
            string property, type, name;
            lineStream >> property >> type >> name;

            string arrayType = TYPE_MAP.count(type) ? TYPE_MAP[type] : "getInt8";
            types[name] = arrayType;
            bytes[name] = SIZE_MAP[type];
            rowOffset += getSizeInBytes(arrayType);
        }
    }

	result.types = types;
	result.bytes = bytes;
	result.byteLength = rowOffset;
	result.numPoints = vertexCount;
	result.headerEndIndex = headerEndIndex;

	auto BoundingBox = getBoundingBox(path,result);
	result.min = BoundingBox.min;
	result.max = BoundingBox.max;

	return result;
}

PlyBox getBoundingBox(string path,PlyHeader header){
	PlyBox result;
	std::ifstream file(path, std::ios::binary);

    if (!file) {
        std::cerr << "Error: Unable to open file " << path << std::endl;
        return result;
    }

    int64_t numRead = 0;
    while (file) {
        // Calculate the current position in the file based on numRead
        int64_t currentPos = header.headerEndIndex + numRead * header.byteLength;

        // Move the file pointer to the current position
        file.seekg(currentPos, std::ios::beg);

        // Buffer to hold 12 bytes for Vector3 (3 floats)
        char buffer[12];

        // Read the first 12 bytes (3 floats) from the current line
        file.read(buffer, 12);

        if (!file) break; // Exit loop if the read fails (e.g., end of file)

        // Convert the buffer to Vector3
		float x, y, z;
		memcpy(&x, &buffer[0], sizeof(float));
		memcpy(&y, &buffer[4], sizeof(float));
		memcpy(&z, &buffer[8], sizeof(float));
		double dx = static_cast<double>(x);
		double dy = static_cast<double>(y);
		double dz = static_cast<double>(z);

        result.min.x = std::min(dx,result.min.x);
		result.min.y = std::min(dy,result.min.y);
		result.min.z = std::min(dz,result.min.z);

		result.max.x = std::max(dx,result.max.x);
		result.max.y = std::max(dy,result.max.y);
		result.max.z = std::max(dz,result.max.z);
		

        // Increment numRead to move to the next position in the file
        numRead++;
    }

    file.close();
	return result;
}