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

constexpr int ATMOS_GRANULARITY = 1;

const float MAX_ATMOS_VELOCITY = std::pow(2.0,ATMOS_GRANULARITY+1)/std::pow(2.0,ATMOS_GRANULARITY);

enum ProcessingLevel {
	TURF_PROCESSING_LEVEL_NONE = -1,
	TURF_PROCESSING_LEVEL_SPACE = 0,
	TURF_PROCESSING_LEVEL_STATION = 1,
	TURF_PROCESSING_LEVEL_PLANET = 2
};

struct ZAtmos {
	ProcessSpecifier processing_turfs;
	std::unordered_map<int,GasMatrix> gasDensities;
	Eigen::SparseMatrix<float> gasEnergy;
	VelocityMatrix gasVelocity;
	Eigen::SparseMatrix<float> energySources;
	VelocityMatrix velocitySources;
	std::unordered_map<int,GasMatrix> sources;
	void resize(int y,int x);
};

class TurfGrid {
public:
	Tile &get(int x, int y, int z);
	Tile &get(int id);
	GasMixture gas_from_turf(int x,int y, int z,bool cache=true);
	GasMixture gas_from_turf(int id,bool cache=true);
	std::complex<float> get_velocity(int id);
	void update_turf(int id, bool can_pass);
	void set_turf_proc_level(int id, int level);
	void add_to_source(int id,GasMixture& gas,float force,float dir);
	void refresh();
	void process();
	void consider_pressure_differences();

private:
	std::vector<ZAtmos> zLevels;
	std::unordered_map<int,GasMixture> gasesByTurf;
	short maxx = 0;
	short maxy = 0;
	short maxz = 0;
	int maxid = 0;
};

GasMixture &get_gas_mixture(Value src);

size_t get_gas_mixture_index(Value val);