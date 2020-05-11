// adapted from monstermos, but not monstermos: it's fastmos

// we're not going quite so far with the atmospherics craziness, instead just moving gas mixes--and ONLY mixes--to C++ land

//i might expand it later with stuff like multithreaded turf reactions, who knows

#include "../core/core.h"
#include "GasMixture.h"
#include "../dmdism/opcodes.h"
#include <chrono>
#include <memory>

trvh fuck(unsigned int args_len, Value* args, Value src)
{
	return Value("fuck");
}

std::unordered_map<std::string, Value> gas_types;
std::unordered_map<unsigned int, int> gas_ids;
//std::unordered_map<unsigned int, std::shared_ptr<GasMixture>> gas_mixtures;
std::vector<Value> gas_id_to_type;
int str_id_extools_pointer;
int gas_mixture_count = 0;
float gas_moles_visible[TOTAL_NUM_GASES];
std::vector<Value> gas_overlays[TOTAL_NUM_GASES];

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
	DecRefCount(args[0].type, args[0].value); // super hacky memory leak fix  (and others) - arguments automatically have refcounts incremented, but not decremented.
	return Value::Null();
}

trvh gasmixture_remove_ratio(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2)
		return Value::Null();
	get_gas_mixture(args[0])->copy_from_mutable(get_gas_mixture(src)->remove_ratio(args[1].valuef));
	DecRefCount(args[0].type, args[0].value);
	return Value::Null();
}

trvh gasmixture_remove(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2)
		return Value::Null();
	get_gas_mixture(args[0])->copy_from_mutable(get_gas_mixture(src)->remove(args[1].valuef));
	DecRefCount(args[0].type, args[0].value);
	return Value::Null();
}

trvh gasmixture_copy_from(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1)
		return Value::Null();
	get_gas_mixture(src)->copy_from_mutable(*get_gas_mixture(args[0]));
	DecRefCount(args[0].type, args[0].value);
	return Value::Null();
}

trvh gasmixture_share(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 1)
		return Value::Null();
	Value ret = Value(get_gas_mixture(src)->share(*get_gas_mixture(args[0]), args_len >= 2 ? args[1].valuef : 4));
	DecRefCount(args[0].type, args[0].value);
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
	get_gas_mixture(src)->set_temperature(args_len > 0 ? args[0].valuef : 0);
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
	return Value::Null();
}

trvh gasmixture_set_moles(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2 || args[0].type != DATUM_TYPEPATH)
		return Value::Null();
	int index = gas_ids[args[0].value];
	get_gas_mixture(src)->set_moles(index, args[1].valuef);
	return Value::Null();
}

trvh gasmixture_scrub_into(unsigned int args_len, Value* args, Value src)
{
	if (args_len < 2)
		return Value::Null();
	// DecRefCount(args[0].type, args[0].value); // Since we're returning it we don't want to decrement the reference count
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
	DecRefCount(args[1].type, args[1].value);
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
void initialize_gas_overlays() {
	Value GLOB = Value::Global().get("GLOB");
	if (!GLOB) return;
	Container meta_gas_info = GLOB.get("meta_gas_info");
	if (!meta_gas_info.type) return;
	for (int i = 0; i < TOTAL_NUM_GASES; ++i)
	{
		Value v = gas_id_to_type[i];
		Container gas_meta = meta_gas_info.at(v);
		gas_moles_visible[i] = gas_meta.at(2);
		gas_overlays[i].clear();
		if (gas_meta.at(3)) {
			Container gas_overlays_list = gas_meta.at(3);
			int num_overlays = gas_overlays_list.length();
			for (int j = 0; j < num_overlays; j++) {
				gas_overlays[i].push_back(gas_overlays_list[j]);
			}
		}
	}
}

int str_id_air;
int str_id_atmosadj;
int str_id_is_openturf;
int str_id_x, str_id_y, str_id_z;
int str_id_current_cycle, str_id_archived_cycle, str_id_planetary_atmos, str_id_initial_gas_mix;
int str_id_active_turfs;
int str_id_react, str_id_consider_pressure_difference, str_id_update_visuals, str_id_floor_rip;

const char* enable_fastmos()
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
	str_id_active_turfs = Core::GetStringId("active_turfs", true);
	str_id_planetary_atmos = Core::GetStringId("planetary_atmos", true);
	str_id_initial_gas_mix = Core::GetStringId("initial_gas_mix", true);
	str_id_react = Core::GetStringId("react", true);
	str_id_consider_pressure_difference = Core::GetStringId("consider pressure difference", true); // byond replaces "_" with " " in proc names. thanks BYOND.
	str_id_update_visuals = Core::GetStringId("update visuals", true);
	str_id_floor_rip = Core::GetStringId("handle decompression floor rip", true);
	str_id_extools_pointer = Core::GetStringId("_extools_pointer_gasmixture", true);

	//Set up gas types map
	std::vector<Value> nullvector = { Value(0.0f) };
	Container gas_types_list = Core::get_proc("/proc/gas_types").call(nullvector);
	Container meta_gas_info = Value::Global().get("meta_gas_info");
	int gaslen = gas_types_list.length();
	if (gaslen != TOTAL_NUM_GASES) {
		return "TOTAL_NUM_GASES does not match the number of /datum/gas subtypes!!";
	}
	for (int i = 0; i < gaslen; ++i)
	{
		Value v = gas_types_list.at(i);
		gas_types[Core::stringify(v)] = gas_types_list.at(i);
		gas_ids[v.value] = i;
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

	return "ok";
}

extern "C" EXPORT const char* init_fastmos(int a, const char** b)
{
	if (!Core::initialize())
	{
		return "Extools Init Failed";
	}
	return enable_fastmos();
}