#pragma once

#include "GasMixture.h"
#include "../core/core.h"
#include <memory>
#include "../third_party/Eigen/Sparse"

class TurfGrid;
struct PlanetAtmosInfo;
struct ExcitedGroup;
struct MonstermosInfo;

extern std::vector<std::weak_ptr<ExcitedGroup>> excited_groups;

struct Tile
{
    Tile();
	void update_air_ref();
	void update_adjacent(TurfGrid &grid);
	void archive(int fire_count);
	void last_share_check();
	void update_planet_atmos();
	// adjacent tiles in the order NORTH,SOUTH,EAST,WEST,UP,DOWN.
	// 
	Tile *adjacent[6];
	unsigned char adjacent_bits = 0;
	size_t air_index;
	Value turf_ref; // not managed because turf refcounts are very unimportant and don't matter
	std::unique_ptr<PlanetAtmosInfo> planet_atmos_info;
	std::shared_ptr<ExcitedGroup> excited_group; // shared_ptr for an actuall good reason this time.
	bool already_processed = false;
};

struct PlanetAtmosInfo
{
	ManagedValue last_initial = Value::Null();
	GasMixture last_mix = GasMixture(monstermos::constants::CELL_VOLUME);
}; // not part of main Tile struct because we don't need it for the whole map

using GasMatrix = Eigen::SparseMatrix<float>;

using VelocityMatrix = Eigen::SparseMatrix<std::complex<float>>;

using ProcessSpecifier = Eigen::Array<char,Eigen::Dynamic,Eigen::Dynamic>;

const int ATMOS_WIDTH = 255;

enum ProcessingLevel {
	PROCESSING_LEVEL_NONE,
	PROCESSING_LEVEL_SPACE,
	PROCESSING_LEVEL_STATION,
	PROCESSING_LEVEL_PLANET
};

class TurfGrid {
public:
	Tile &get(int x, int y, int z);
	Tile &get(int id);
	GasMixture gas_from_turf(int x,int y, int z,bool cache=true);
	GasMixture gas_from_turf(int id,bool cache=true);
	void refresh();
	void process();

private:
	std::vector<ProcessSpecifier> processing_turfs;
	std::unordered_map<int,std::vector<GasMatrix>> gasDensities;
	std::vector<GasMatrix> gasEnergy;
	std::vector<VelocityMatrix> gasVelocity;
	std::vector<GasMatrix> energySources;
	std::vector<VelocityMatrix> velocitySources;
	std::unordered_map<int,std::vector<GasMatrix>> sources;
	std::unordered_map<int,GasMixture> gasesByTurf;
	short maxx = 0;
	short maxy = 0;
	short maxz = 0;
	int maxid = 0;
};

GasMixture &get_gas_mixture(Value src);

size_t get_gas_mixture_index(Value val);