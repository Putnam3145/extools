#include "monstermos.h"

#include "../core/core.h"
#include "GasMixture.h"
#include "Reaction.h"
#include "turf_grid.h"
#include "../dmdism/opcodes.h"

#include <cmath>
#include <chrono>
#include <algorithm>
#include <list>
#include <unordered_set>
#include <mutex>
#include <utility>
#include <tuple>
#include <thread>

using namespace monstermos::constants;

trvh fuck(unsigned int args_len, Value* args, Value src)
{
	return Value("fuck");
}

std::unordered_map<std::string, Value> gas_types;
std::unordered_map<unsigned int, int> gas_ids;
//std::unordered_map<unsigned int, std::shared_ptr<GasMixture>> gas_mixtures;
std::vector<Value> gas_id_to_type;
std::vector<std::shared_ptr<Reaction>> cached_reactions;
TurfGrid all_turfs;
Value SSair;
int str_id_extools_pointer;
int gas_mixture_count = 0;
float gas_moles_visible[TOTAL_NUM_GASES];
float gas_fusion_power[TOTAL_NUM_GASES];
std::vector<Value> gas_overlays[TOTAL_NUM_GASES];

int o2,plasma,tritium,co2,water_vapor,n2o,bz,no2;

std::shared_ptr<GasMixture> &get_gas_mixture(Value val)
{
	uint32_t v = val.get_by_id(str_id_extools_pointer).value;
	if (v == 0) Runtime("Gas mixture has null extools pointer");
	return *((std::shared_ptr<GasMixture>*)v);
}

int str_id_volume;
trvh gasmixture_register(unsigned int args_len, Value* args, Value src)
{
	//gas_mixtures[src.value] = std::make_shared<GasMixture>(src.get_by_id(str_id_volume).valuef);
	std::shared_ptr<GasMixture> *ptr = new std::shared_ptr<GasMixture>;
	*ptr = std::make_shared<GasMixture>(src.get_by_id(str_id_volume).valuef);
	SetVariable(src.type, src.value, str_id_extools_pointer, Value(NUMBER, (int)ptr));
	gas_mixture_count++;
	return Value::Null();
}

trvh gasmixture_unregister(unsigned int args_len, Value* args, Value src)
{
	uint32_t v = src.get_by_id(str_id_extools_pointer).value;
	if (v != 0) {
		std::shared_ptr<GasMixture> *gm = (std::shared_ptr<GasMixture> *)v;
		delete gm;
		gas_mixture_count--;
		SetVariable(src.type, src.value, str_id_extools_pointer, Value::Null());
	}
	return Value::Null();
}

DelDatumPtr oDelDatum;
void hDelDatum(unsigned int datum_id) {
	RawDatum *datum = Core::GetDatumPointerById(datum_id);
	if (datum != nullptr) {
		std::shared_ptr<GasMixture> *gm = nullptr;
		if (datum->len_vars < 10) { // if it has a whole bunch of vars it's probably not a gas mixture. Please don't add a whole bunch of vars to gas mixtures.
			for (int i = 0; i < datum->len_vars; i++) {
				if (datum->vars[i].id == str_id_extools_pointer) {
					gm = (std::shared_ptr<GasMixture> *)datum->vars[i].value.value;
					datum->vars[i].value = Value::Null();
					break;
				}
			}
		}
		if (gm != nullptr) {
			delete gm;
			gas_mixture_count--;
		}
	}
	oDelDatum(datum_id);
}

trvh gasmixture_heat_capacity(unsigned int args_len, Value* args, Value src)
{
	return Value(get_gas_mixture(src)->heat_capacity());
}

trvh gasmixture_set_min_heat_capacity(unsigned int args_len, Value* args, Value src)
{
	get_gas_mixture(src)->set_min_heat_capacity(args_len > 0 ? args[0].valuef : 0);
	return Value::Null();
}

trvh gasmixture_total_moles(unsigned int args_len, Value* args, Value src)
{
	return Value(get_gas_mixture(src)->total_moles());
}

trvh gasmixture_return_pressure(unsigned int args_len, Value* args, Value src)
{
	return Value(get_gas_mixture(src)->return_pressure());
}

trvh gasmixture_return_temperature(unsigned int args_len, Value* args, Value src)
{
	return Value(get_gas_mixture(src)->get_temperature());
}

trvh gasmixture_return_volume(unsigned int args_len, Value* args, Value src)
{
	return Value(get_gas_mixture(src)->get_volume());
}

trvh gasmixture_thermal_energy(unsigned int args_len, Value* args, Value src)
{
	return Value(get_gas_mixture(src)->thermal_energy());
}

trvh gasmixture_archive(unsigned int args_len, Value* args, Value src)
{
	get_gas_mixture(src)->archive();
	return Value::Null();
}

trvh gasmixture_merge(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1)
		return Value::Null();
	get_gas_mixture(src)->merge(*get_gas_mixture(args[0]));
	return Value::Null();
}

trvh gasmixture_remove_ratio(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2)
		return Value::Null();
	get_gas_mixture(args[0])->copy_from_mutable(get_gas_mixture(src)->remove_ratio(args[1].valuef));
	return Value::Null();
}

trvh gasmixture_remove(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2)
		return Value::Null();
	get_gas_mixture(args[0])->copy_from_mutable(get_gas_mixture(src)->remove(args[1].valuef));
	return Value::Null();
}

trvh gasmixture_copy_from(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1)
		return Value::Null();
	get_gas_mixture(src)->copy_from_mutable(*get_gas_mixture(args[0]));
	return Value::Null();
}

trvh gasmixture_share(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1)
		return Value::Null();
	Value ret = Value(get_gas_mixture(src)->share(*get_gas_mixture(args[0]), args_len >= 2 ? args[1].valuef : 4));
	return ret;
}

trvh gasmixture_get_last_share(unsigned int args_len, Value* args, Value src)
{
	return Value(get_gas_mixture(src)->get_last_share());
}

trvh gasmixture_get_gases(unsigned int args_len, Value* args, Value src)
{
	List l(CreateList(0));
	GasMixture &gm = *get_gas_mixture(src);
	for (int i = 0; i < TOTAL_NUM_GASES; i++) {
		if (gm.get_moles(i) >= GAS_MIN_MOLES) {
			l.append(gas_id_to_type[i]);
		}
	}
	return l;
}

trvh gasmixture_set_temperature(unsigned int args_len, Value* args, Value src)
{
	float vf = args_len > 0 ? args[0].valuef : 0;
	if (std::isnan(vf) || std::isinf(vf)) {
		Runtime("Attempt to set temperature to NaN or Infinity");
	} else {
		GasMixture &src_gas = *get_gas_mixture(src);
		src_gas.set_temperature(vf);
		src_gas.set_dirty(true);
	}
	return Value::Null();
}

trvh gasmixture_set_volume(unsigned int args_len, Value* args, Value src)
{
	get_gas_mixture(src)->set_volume(args_len > 0 ? args[0].valuef : 0);
	return Value::Null();
}

trvh gasmixture_get_moles(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1 || args[0].type != DATUM_TYPEPATH)
		return Value::Null();
	int index = gas_ids[args[0].value];
	return Value(get_gas_mixture(src)->get_moles(index));
}

trvh gasmixture_set_moles(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2 || args[0].type != DATUM_TYPEPATH)
		return Value::Null();
	int index = gas_ids[args[0].value];
	GasMixture &src_gas = *get_gas_mixture(src);
	src_gas.set_moles(index, args[1].valuef);
	src_gas.set_dirty(true);
	return Value::Null();
}

trvh gasmixture_scrub_into(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2)
		return Value::Null();
	GasMixture &src_gas = *get_gas_mixture(src);
	GasMixture &dest_gas = *get_gas_mixture(args[0]);
	Container gases_to_scrub = args[1];
	int num_gases = gases_to_scrub.length();
	GasMixture buffer(CELL_VOLUME);
	buffer.set_temperature(src_gas.get_temperature());
	for (int i = 0; i < num_gases; i++) {
		Value typepath = gases_to_scrub[i];
		if (typepath.type != DATUM_TYPEPATH) continue;
		int index = gas_ids[typepath.value];
		buffer.set_moles(index, buffer.get_moles(index) + src_gas.get_moles(index));
		src_gas.set_moles(index, 0);
	}
	dest_gas.merge(buffer);
	IncRefCount(args[0].type, args[0].value);
	return args[0];
}

trvh gasmixture_mark_immutable(unsigned int args_len, Value* args, Value src)
{
	get_gas_mixture(src)->mark_immutable();
	return Value::Null();
}

trvh gasmixture_clear(unsigned int args_len, Value* args, Value src)
{
	get_gas_mixture(src)->clear();
	return Value::Null();
}

trvh gasmixture_compare(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1)
		return Value::Null();
	int result = get_gas_mixture(src)->compare(*get_gas_mixture(args[0]));
	if (result == -1) {
		return Value("temp");
	}
	else if (result == -2) {
		return Value("");
	} else{
		return gas_id_to_type[result];
	}
}

trvh gasmixture_multiply(unsigned int args_len, Value* args, Value src)
{
	get_gas_mixture(src)->multiply(args_len > 0 ? args[0].valuef : 1);
	return Value::Null();
}

trvh gasmixture_react(unsigned int args_len, Value* args, Value src)
{
	GasMixture &src_gas = *get_gas_mixture(src);
	if(src_gas.sleeping()) return Value((float)NO_REACTION);
	auto ret = 0;
	Value holder;
	if(args_len == 0)
	{
		holder = Value();
	}
	else
	{
		holder = args[0];
	}
	for(int i=0;i<cached_reactions.size();i++)
	{
		auto reaction = cached_reactions[i];
		if(reaction->check_conditions(src_gas)) {
			IncRefCount(src.type,src.value); // have to do this or the gas mixture will be GC'd at the end of the function
			IncRefCount(holder.type,holder.value); // i'm assuming this would also end up GC'd--even worse
			ret |= cached_reactions[i]->react(src_gas,src,holder);
		}
		if(ret & STOP_REACTIONS) return Value((float)ret);
	}
	src_gas.set_dirty(ret != NO_REACTION);
	return Value((float)ret);
}

trvh turf_update_adjacent(unsigned int args_len, Value* args, Value src)
{
	if (src.type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(src.value);
	if (tile != nullptr) {
		tile->update_adjacent(all_turfs);
	}
	return Value::Null();
}

trvh turf_update_air_ref(unsigned int args_len, Value* args, Value src)
{
	if (src.type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(src.value);
	if (tile != nullptr) {
		tile->update_air_ref();
	}
	return Value::Null();
}

trvh turf_eg_reset_cooldowns(unsigned int args_len, Value* args, Value src)
{
	if (src.type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(src.value);
	if (tile != nullptr) {
		if (tile->excited_group) {
			tile->excited_group->reset_cooldowns();
		}
	}
	return Value::Null();
}

trvh turf_eg_garbage_collect(unsigned int args_len, Value* args, Value src)
{
	if (src.type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(src.value);
	if (tile != nullptr) {
		if (tile->excited_group) {
			// store to local variable to prevent it from being destructed while we're still using it because that causes segfaults.
			std::shared_ptr<ExcitedGroup> eg = tile->excited_group;
			eg->dismantle(false);
		}
	}
	return Value::Null();
}

trvh turf_get_excited(unsigned int args_len, Value* args, Value src)
{
	if (src.type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(src.value);
	if (tile != nullptr) {
		return Value(tile->excited ? 1.0 : 0.0);
	}
	return Value::Null();
}
trvh turf_set_excited(unsigned int args_len, Value* args, Value src)
{
	if (src.type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(src.value);
	if (tile != nullptr) {
		tile->excited = args_len > 0 ? (bool)args[0] : false;
	}
	return Value::Null();
}

trvh turf_eq(unsigned int args_len, Value* args, Value src) {
	if (src.type != TURF || args_len < 1) { return Value::Null(); }
	Tile *tile = all_turfs.get(src.value);
	if (tile != nullptr) {
		tile->equalize_pressure_in_zone(args[0]);
	}
	return Value::Null();
}

unsigned int str_id_atmos_overlay_types;
trvh turf_update_visuals(unsigned int args_len, Value* args, Value src) {
	if (src.type != TURF) { return Value::Null(); }
	Tile* tile = all_turfs.get(src.value);
	if (!tile->air) return Value::Null();
	GasMixture& gm = *tile->air;
	Value old_overlay_types_val = src.get_by_id(str_id_atmos_overlay_types);
	std::vector<Value> overlay_types;

	for (int i = 0; i < TOTAL_NUM_GASES; i++) {
		if (!gas_overlays[i].size()) continue;
		if (gm.get_moles(i) > gas_moles_visible[i]) {
			// you know whats fun?
			// getting cucked by BYOND arrays starting at 1. How did this not segfault before? Beats me! I love undefined behavior!    Bandaid: VV
			overlay_types.push_back(gas_overlays[i][std::fmin(FACTOR_GAS_VISIBLE_MAX, (int)std::ceil(gm.get_moles(i) / MOLES_GAS_VISIBLE_STEP))-1]);
		}
	}

	if (!overlay_types.size() && !old_overlay_types_val) return Value::Null();
	if (old_overlay_types_val) {
		List old_overlay_types(old_overlay_types_val);
		if (overlay_types.size() == old_overlay_types.list->length) {
			bool is_different = false;
			for (int i = 0; i < overlay_types.size(); i++) {
				if (overlay_types[i] != old_overlay_types.at(i)) {
					is_different = true; break;
				}
			}
			if (!is_different) {
				return Value::Null();
			}
		}
	}
	
	List l(CreateList(0));
	for (int i = 0; i < overlay_types.size(); i++) {
		l.append(overlay_types[i]);
	}
	src.invoke("set_visuals", { Value(l) } );
	return Value::Null();
}

class Stopwatch {
	private:
		std::chrono::time_point<std::chrono::steady_clock> start;
	public:
		Stopwatch() {
			start = std::chrono::steady_clock::now();
		}
		void restart() {
			start = std::chrono::steady_clock::now();
		}
		int peek() {
			auto end = std::chrono::steady_clock::now();
			return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
		}
};

std::vector<std::weak_ptr<ExcitedGroup>> excited_groups_currentrun;
trvh SSair_process_excited_groups(unsigned int args_len, Value* args, Value src) {
	auto sw = Stopwatch();
	float time_limit = args[1].valuef * 100000.0f;

	if (args_len < 2) { return Value::Null(); }
	if (!args[0]) {
		excited_groups_currentrun = excited_groups; // this copies it.... right?
	}
	while (excited_groups_currentrun.size()) {
		std::shared_ptr<ExcitedGroup> eg = excited_groups_currentrun.back().lock();
		excited_groups_currentrun.pop_back();
		if (!eg) continue;
		eg->breakdown_cooldown++;
		eg->dismantle_cooldown++;
		if (eg->breakdown_cooldown >= EXCITED_GROUP_BREAKDOWN_CYCLES)
			eg->self_breakdown();
		if (eg->dismantle_cooldown >= EXCITED_GROUP_DISMANTLE_CYCLES)
			eg->dismantle(true);
		if (sw.peek() > time_limit) {
			return Value::True();
		}
	}
	return Value::False();
}

trvh SSair_get_amt_excited_groups(unsigned int args_len, Value* args, Value src) {
	return Value(excited_groups.size());
}

std::mutex process_mutex;
std::mutex done_mutex;
std::unordered_set<Tile*> active_turfs;
std::list<Tile*> active_turfs_currentrun;
std::list< std::tuple< Tile*, std::vector<Tile*>, int, bool > > processing_turfs;
std::list<std::pair < Tile*, std::vector< std::pair<Tile*, float> > > > done_processing_turfs;

void add_to_active(Tile* t)
{
	active_turfs.insert(t);
}

void remove_from_active(Tile* t)
{
	active_turfs.erase(t);
}

trvh SSair_add_to_active(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1 || args[0].type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(args[0].value);
	if (tile != nullptr) {
		add_to_active(tile);
	}
	return Value::Null();
}

trvh SSair_remove_from_active(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1 || args[0].type != TURF) { return Value::Null(); }
	Tile *tile = all_turfs.get(args[0].value);
	if (tile != nullptr) {
		remove_from_active(tile);
	}
	return Value::Null();
}

trvh SSair_clear_active_turfs(unsigned int args_len, Value* args, Value src)
{
	active_turfs.clear();
	return Value::Null();
}

trvh SSair_active_turf_length(unsigned int args_len, Value* args, Value src)
{
	return Value((float)(active_turfs.size()));
}

trvh SSair_active_run_length(unsigned int args_len, Value* args, Value src)
{
	return Value((float)(active_turfs_currentrun.size()));
}

trvh SSair_processing_length(unsigned int args_len, Value* args, Value src)
{
	return Value((float)(processing_turfs.size()));
}

trvh SSair_done_length(unsigned int args_len, Value* args, Value src)
{
	return Value((float)(done_processing_turfs.size()));
}

#include <algorithm>

trvh SSair_get_active_turfs(unsigned int args_len, Value* args, Value src)
{
	List l(CreateList(0));
	std::for_each(active_turfs.begin(),active_turfs.end(),[&](Tile* t) mutable {
		l.append(t->turf_ref);
	});
	return l;
}

trvh ssair_process_active_turfs(unsigned int args_len, Value* args, Value src)
{
	auto sw = Stopwatch();
	float time_limit = args[1].valuef * 100000.0f;
	int fire_count = SSair.get("times_fired");
	if(active_turfs_currentrun.size() == 0) {
		std::for_each(active_turfs.begin(),active_turfs.end(),[](Tile* t) {
			active_turfs_currentrun.push_back(t);
		});
	}
	while(active_turfs_currentrun.size()) {
		Tile& cur_turf = *active_turfs_currentrun.front();
		active_turfs_currentrun.pop_front();
		auto res = cur_turf.pre_process_cell(fire_count);
		if(std::get<1>(res) >= 0)
		{
			std::unique_lock<std::mutex> lock(process_mutex);
			processing_turfs.push_back(std::tuple_cat(std::make_tuple(&cur_turf),res));
		}
		if (sw.peek() > time_limit) {
			return Value::True();
		}
	}
	return Value::False();
}

bool continue_processing_atmos = true;

void air_process_loop()
{
	while(continue_processing_atmos)
	{
		if(processing_turfs.size())
		{
			process_mutex.lock();
			auto payload = processing_turfs.front();
			processing_turfs.pop_front();
			process_mutex.unlock();
			auto [ tile, enemy_tiles, adjacent_tiles, has_planetary_atmos ] = payload;
			auto differences = tile->process_cell(enemy_tiles,adjacent_tiles,has_planetary_atmos);
			{
				std::unique_lock<std::mutex> lock(done_mutex);
				done_processing_turfs.push_back(make_pair(tile,differences));
			}
		}
		else
		{
			std::this_thread::yield();
		}
	}
}

trvh ssair_post_process_turfs(unsigned int args_len, Value* args, Value src)
{
	auto sw = Stopwatch();
	float time_limit = args[1] * 100000.0f;
	while(done_processing_turfs.size()) {
		done_mutex.lock();
		auto cur_turf = done_processing_turfs.front();
		done_processing_turfs.pop_front();
		done_mutex.unlock();
		cur_turf.first->post_process_cell(cur_turf.second);
		if (sw.peek() > time_limit) {
			return Value::True();
		}
	}
	return Value::False();
}

trvh ssair_process_turf_equalize(unsigned int args_len, Value* args, Value src)
{
	auto sw = Stopwatch();
	float time_limit = args[1].valuef * 100000.0f;
	int fire_count = SSair.get("times_fired");
	if(active_turfs_currentrun.size() == 0) {
		std::for_each(active_turfs.begin(),active_turfs.end(),[](Tile* t) {
			active_turfs_currentrun.push_back(t);
		});
	}
	while(active_turfs_currentrun.size()) {
		auto cur_turf = active_turfs_currentrun.front();
		active_turfs_currentrun.pop_front();
		cur_turf->equalize_pressure_in_zone(fire_count);
		if (sw.peek() > time_limit) {
			return Value::True();
		}
	}
	return Value::False();
}

trvh refresh_atmos_grid(unsigned int args_len, Value* args, Value src)
{
	all_turfs.refresh();
	return Value::Null();
}

void initialize_gas_overlays() {
	Value GLOB = Value::Global().get("GLOB");
	if (!GLOB) return;
	Container meta_gas_visibility = GLOB.get("meta_gas_visibility");
	Container meta_gas_overlays = GLOB.get("meta_gas_overlays");
	Container meta_gas_fusions = GLOB.get("meta_gas_fusions");
	if (!meta_gas_visibility.type) return;
	for (int i = 0; i < TOTAL_NUM_GASES; ++i)
	{
		Value v = gas_id_to_type[i];
		gas_moles_visible[i] = meta_gas_visibility.at(v);
		gas_overlays[i].clear();
		gas_fusion_power[i] = meta_gas_fusions.at(v);
		Container gas_overlays_list = meta_gas_overlays.at(v);
		int num_overlays = gas_overlays_list.length();
		for (int j = 0; j < num_overlays; j++) {
			gas_overlays[i].push_back(gas_overlays_list[j]);
		}
	}
}

trvh SSair_update_ssair(unsigned int args_len, Value* args, Value src) {
	SSair = src;
	initialize_gas_overlays();
	return Value::Null();
}

trvh SSair_update_gas_reactions(unsigned int args_len, Value* args, Value src) {
	Container gas_reactions = SSair.get("gas_reactions");
	cached_reactions.clear();
	for(int i = 0; i < gas_reactions.length(); i++)
	{
		cached_reactions.push_back(std::make_shared<ByondReaction>(gas_reactions.at(i)));
	}
	std::sort(cached_reactions.begin(),cached_reactions.end(),
	[](std::shared_ptr<Reaction> a, std::shared_ptr<Reaction> b) { return a->get_priority() > b->get_priority(); });
	return Value::Null();
}

int str_id_air;
int str_id_atmosadj;
int str_id_is_openturf;
int str_id_x, str_id_y, str_id_z;
int str_id_current_cycle, str_id_archived_cycle, str_id_planetary_atmos, str_id_initial_gas_mix;
int str_id_react, str_id_gas_reactions, str_id_consider_pressure_difference, str_id_update_visuals, str_id_floor_rip;
int str_id_monstermos_turf_limit, str_id_monstermos_hard_turf_limit;

std::thread processing_loop;

trvh end_thread(unsigned int args_len, Value* args, Value src) {
	continue_processing_atmos = false;
	processing_loop.join();
	return Value::Null();
}

const char* enable_monstermos()
{
	oDelDatum = (DelDatumPtr)Core::install_hook((void*)DelDatum, (void*)hDelDatum);
	// get the var IDs for SANIC SPEED
	str_id_air = Core::GetStringId("air", true);
	str_id_atmosadj = Core::GetStringId("atmos_adjacent_turfs", true);
	str_id_volume = Core::GetStringId("initial_volume", true);
	str_id_is_openturf = Core::GetStringId("is_openturf", true);
	str_id_x = Core::GetStringId("x", true);
	str_id_y = Core::GetStringId("y", true);
	str_id_z = Core::GetStringId("z", true);
	str_id_current_cycle = Core::GetStringId("current_cycle", true);
	str_id_archived_cycle = Core::GetStringId("archived_cycle", true);
	str_id_planetary_atmos = Core::GetStringId("planetary_atmos", true);
	str_id_initial_gas_mix = Core::GetStringId("initial_gas_mix", true);
	str_id_atmos_overlay_types = Core::GetStringId("atmos_overlay_types", true);
	str_id_react = Core::GetStringId("react", true);
	str_id_gas_reactions = Core::GetStringId("gas_reactions", true);
	str_id_consider_pressure_difference = Core::GetStringId("consider pressure difference", true); // byond replaces "_" with " " in proc names. thanks BYOND.
	str_id_update_visuals = Core::GetStringId("update visuals", true);
	str_id_floor_rip = Core::GetStringId("handle decompression floor rip", true);
	str_id_extools_pointer = Core::GetStringId("_extools_pointer_gasmixture", true);
	str_id_monstermos_turf_limit = Core::GetStringId("monstermos_turf_limit",true);
	str_id_monstermos_hard_turf_limit = Core::GetStringId("monstermos_turf_limit",true);

	SSair = Value::Global().get("SSair");
	//Set up gas types map
	std::vector<Value> nullvector = { Value(0.0f) };
	Container gas_types_list = Core::get_proc("/proc/gas_types").call(nullvector);
	int gaslen = gas_types_list.length();
	if (gaslen != TOTAL_NUM_GASES) {
		return "TOTAL_NUM_GASES does not match the number of /datum/gas subtypes!!";
	}
	for (int i = 0; i < gaslen; ++i)
	{
		Value v = gas_types_list.at(i);
		std::string type_name = Core::stringify(v);
		gas_types[type_name] = gas_types_list.at(i);
		gas_ids[v.value] = i;
		if(type_name == "/datum/gas/oxygen") o2 = i;
		else if(type_name == "/datum/gas/carbon_dioxide") co2 = i;
		else if(type_name == "/datum/gas/tritium") tritium = i;
		else if(type_name == "/datum/gas/plasma") plasma = i;
		else if(type_name == "/datum/gas/water_vapor") water_vapor = i;
		else if(type_name == "/datum/gas/nitrous_oxide") n2o = i;
		else if(type_name == "/datum/gas/nitryl") no2 = i;
		else if(type_name == "/datum/gas/bz") bz = i;
		gas_specific_heat[i] = gas_types_list.at(v).valuef;
		gas_id_to_type.push_back(v);
	}
	initialize_gas_overlays();
	//Set up hooks
	Core::get_proc("/datum/gas_mixture/proc/__gasmixture_register").hook(gasmixture_register);
	Core::get_proc("/datum/gas_mixture/proc/__gasmixture_unregister").hook(gasmixture_unregister);
	Core::get_proc("/datum/gas_mixture/proc/heat_capacity").hook(gasmixture_heat_capacity);
	Core::get_proc("/datum/gas_mixture/proc/set_min_heat_capacity").hook(gasmixture_set_min_heat_capacity);
	Core::get_proc("/datum/gas_mixture/proc/total_moles").hook(gasmixture_total_moles);
	Core::get_proc("/datum/gas_mixture/proc/return_pressure").hook(gasmixture_return_pressure);
	Core::get_proc("/datum/gas_mixture/proc/return_temperature").hook(gasmixture_return_temperature);
	Core::get_proc("/datum/gas_mixture/proc/return_volume").hook(gasmixture_return_volume);
	Core::get_proc("/datum/gas_mixture/proc/thermal_energy").hook(gasmixture_thermal_energy);
	Core::get_proc("/datum/gas_mixture/proc/archive").hook(gasmixture_archive);
	Core::get_proc("/datum/gas_mixture/proc/merge").hook(gasmixture_merge);
	Core::get_proc("/datum/gas_mixture/proc/copy_from").hook(gasmixture_copy_from);
	Core::get_proc("/datum/gas_mixture/proc/share").hook(gasmixture_share);
	Core::get_proc("/datum/gas_mixture/proc/compare").hook(gasmixture_compare);
	Core::get_proc("/datum/gas_mixture/proc/get_gases").hook(gasmixture_get_gases);
	Core::get_proc("/datum/gas_mixture/proc/__remove").hook(gasmixture_remove);
	Core::get_proc("/datum/gas_mixture/proc/__remove_ratio").hook(gasmixture_remove_ratio);
	Core::get_proc("/datum/gas_mixture/proc/set_temperature").hook(gasmixture_set_temperature);
	Core::get_proc("/datum/gas_mixture/proc/set_volume").hook(gasmixture_set_volume);
	Core::get_proc("/datum/gas_mixture/proc/get_moles").hook(gasmixture_get_moles);
	Core::get_proc("/datum/gas_mixture/proc/set_moles").hook(gasmixture_set_moles);
	Core::get_proc("/datum/gas_mixture/proc/scrub_into").hook(gasmixture_scrub_into);
	Core::get_proc("/datum/gas_mixture/proc/mark_immutable").hook(gasmixture_mark_immutable);
	Core::get_proc("/datum/gas_mixture/proc/clear").hook(gasmixture_clear);
	Core::get_proc("/datum/gas_mixture/proc/multiply").hook(gasmixture_multiply);
	Core::get_proc("/datum/gas_mixture/proc/get_last_share").hook(gasmixture_get_last_share);
	Core::get_proc("/datum/gas_mixture/proc/react").hook(gasmixture_react);
	Core::get_proc("/turf/proc/__update_extools_adjacent_turfs").hook(turf_update_adjacent);
	Core::get_proc("/turf/proc/update_air_ref").hook(turf_update_air_ref);
	Core::get_proc("/turf/open/proc/eg_reset_cooldowns").hook(turf_eg_reset_cooldowns);
	Core::get_proc("/turf/open/proc/eg_garbage_collect").hook(turf_eg_garbage_collect);
	Core::get_proc("/turf/open/proc/get_excited").hook(turf_get_excited);
	Core::get_proc("/turf/open/proc/set_excited").hook(turf_set_excited);
	Core::get_proc("/turf/open/proc/equalize_pressure_in_zone").hook(turf_eq);
	Core::get_proc("/turf/open/proc/update_visuals").hook(turf_update_visuals);
	Core::get_proc("/world/proc/refresh_atmos_grid").hook(refresh_atmos_grid);
	Core::get_proc("/datum/controller/subsystem/air/proc/process_excited_groups_extools").hook(SSair_process_excited_groups);
	Core::get_proc("/datum/controller/subsystem/air/proc/process_turf_equalize_extools").hook(ssair_process_turf_equalize);
	Core::get_proc("/datum/controller/subsystem/air/proc/process_active_turfs_extools").hook(ssair_process_active_turfs);
	Core::get_proc("/datum/controller/subsystem/air/proc/post_process_turfs_extools").hook(ssair_post_process_turfs);
	Core::get_proc("/datum/controller/subsystem/air/proc/get_amt_excited_groups").hook(SSair_get_amt_excited_groups);
	Core::get_proc("/datum/controller/subsystem/air/proc/extools_update_ssair").hook(SSair_update_ssair);
	Core::get_proc("/datum/controller/subsystem/air/proc/extools_setup_gas_reactions").hook(SSair_update_gas_reactions);
	Core::get_proc("/datum/controller/subsystem/air/proc/remove_from_active_extools").hook(SSair_remove_from_active);
	Core::get_proc("/datum/controller/subsystem/air/proc/add_to_active_extools").hook(SSair_add_to_active);
	Core::get_proc("/datum/controller/subsystem/air/proc/clear_active_turfs").hook(SSair_clear_active_turfs);
	Core::get_proc("/datum/controller/subsystem/air/proc/active_turfs_length").hook(SSair_active_turf_length);
	Core::get_proc("/datum/controller/subsystem/air/proc/cpp_currentrun_length").hook(SSair_active_run_length);
	Core::get_proc("/datum/controller/subsystem/air/proc/processing_length").hook(SSair_processing_length);
	Core::get_proc("/datum/controller/subsystem/air/proc/post_processing_length").hook(SSair_done_length);
	Core::get_proc("/datum/controller/subsystem/air/proc/get_active_turfs").hook(SSair_get_active_turfs);
	Core::get_proc("/proc/destroy_extools_atmos_thread").hook(end_thread);
	processing_loop = std::thread(air_process_loop);
	all_turfs.refresh();
	return "ok";
}
