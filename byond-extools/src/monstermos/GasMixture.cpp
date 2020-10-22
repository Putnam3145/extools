#include "GasMixture.h"

#include <algorithm>
#include <cstring>
#include <cmath>

#include <vector>

#include <numeric>

#include <execution>

using namespace monstermos::constants;

std::vector<float> gas_specific_heat;

extern int total_num_gases;

void GasMixture::garbage_collect() {
    int last_non_zero = 0;
    for(int i=0;i<moles.size();i++)
    {
        if(moles[i] >= GAS_MIN_MOLES)
        {
            last_non_zero = i;
        }
    }
    moles.resize(last_non_zero+1);
}

GasMixture::GasMixture(float v)
{
    if(v < 0) v = 0;
    volume = v;
    moles.resize(total_num_gases);
    moles_archived.resize(total_num_gases);
}

GasMixture::GasMixture() {}

float GasMixture::get_moles(int gas_type) const {
    return moles.size() > gas_type ? moles[gas_type] : 0.0;
}

void GasMixture::set_moles(int gas_type, float new_moles) { 
    if(!immutable) {
        if(gas_type > moles.size())
        {
            moles.resize(gas_type+1);
        }
        moles[gas_type] = new_moles;
    }
}

float GasMixture::total_moles() const {
    return std::reduce(std::execution::seq,moles.cbegin(),moles.cend());
}

void GasMixture::mark_immutable() {
    immutable = true;
}

float GasMixture::heat_capacity() const {
    return std::max(
        (float)std::transform_reduce(
            std::execution::seq,
            moles.cbegin(),moles.cend(),
            gas_specific_heat.cbegin(),
            0.0),
        min_heat_capacity);
}

float GasMixture::heat_capacity_archived() const {
    return std::max(
        (float)std::transform_reduce(
            std::execution::seq,
            moles_archived.cbegin(),moles_archived.cend(),
            gas_specific_heat.cbegin(),
            0.0), 
        min_heat_capacity);
}

void GasMixture::set_min_heat_capacity(float n) {
	if (immutable) return;
	min_heat_capacity = n;
}

float GasMixture::return_pressure() const {
    if(volume <= 0) return 0;
    return total_moles() * R_IDEAL_GAS_EQUATION * temperature / volume;
}

float GasMixture::thermal_energy() const {
    return temperature * heat_capacity();
}

void GasMixture::archive() {
    moles_archived = moles;
    temperature_archived = temperature;
}

void GasMixture::merge(const GasMixture &giver) {
    if(immutable) return;
    if(std::abs(temperature - giver.temperature) > MINIMUM_TEMPERATURE_DELTA_TO_CONSIDER) {
        float self_heat_capacity = heat_capacity();
        float giver_heat_capacity = giver.heat_capacity();
        float combined_heat_capacity = self_heat_capacity + giver_heat_capacity;
        if(combined_heat_capacity) {
			temperature = (giver.temperature * giver_heat_capacity + temperature * self_heat_capacity) / combined_heat_capacity;
        }
    }
    moles.resize(std::max(moles.size(),giver.moles.size()));
    std::transform(
        std::execution::seq,
        moles.begin(),moles.end(),
        giver.moles.cbegin(),
        moles.begin(),
        std::plus<float>());
}

GasMixture GasMixture::remove(float amount) {
	return remove_ratio(amount / total_moles());
}

GasMixture GasMixture::remove_ratio(float ratio) {
    if(ratio <= 0)
        return GasMixture(volume);

    if(ratio > 1) ratio = 1;

    auto removed = GasMixture(volume);
    removed.temperature = temperature;
    std::transform(std::execution::seq,
        moles.begin(),moles.end(),
        removed.moles.begin(),
        [&ratio](auto& gas) {
            return gas * ratio;
    });
    if(!immutable) [[likely]]
    {
        std::transform(std::execution::seq,
        moles.begin(),moles.end(),
        removed.moles.begin(),
        moles.begin(),
        [&ratio](auto& myGas,auto& removedGas) {
            return myGas-removedGas;
        });
    }
    removed.garbage_collect();
    garbage_collect();
    return removed;
}

void GasMixture::copy_from_mutable(const GasMixture &sample) {
    if(immutable) return;
    moles = sample.moles;
    temperature = sample.temperature;
}

float GasMixture::share(GasMixture &sharer, const int atmos_adjacent_turfs) {
    const float temperature_delta = temperature_archived - sharer.temperature_archived;
    const float abs_temperature_delta = std::abs(temperature_delta);
    const float old_self_heat_capacity = heat_capacity();
    const float old_sharer_heat_capacity = sharer.heat_capacity();
    float heat_capacity_self_to_sharer = 0;
    float heat_capacity_sharer_to_self = 0;
    float moved_moles = 0;
    float abs_moved_moles = 0;
    auto capacities = gas_specific_heat; // cache for sanic speed (yeah, seriously)
    const int max_size = std::max(std::max(moles_archived.size(),sharer.moles_archived.size()),std::max(moles.size(),sharer.moles.size()));
    moles.resize(max_size);
    sharer.moles.resize(max_size);
    moles_archived.resize(max_size);
    sharer.moles_archived.resize(max_size);
    auto moles_copy = moles;
    auto sharer_copy = sharer.moles;
    std::vector<float> heat_cap_deltas(max_size);
    for(int i = 0; i < max_size; i++) {
        const float delta = (moles_archived[i] - sharer.moles_archived[i])/(atmos_adjacent_turfs+1);
        heat_cap_deltas[i] = delta * capacities[i];
        moles_copy[i] -= delta;
        sharer_copy[i] += delta;
        moved_moles += delta;
        abs_moved_moles += std::abs(delta);
    }
    if((abs_temperature_delta > MINIMUM_TEMPERATURE_DELTA_TO_CONSIDER)) {
        for(int i=0;i<max_size;i++) {
            const float gas_heat_capacity = heat_cap_deltas[i];
            if(gas_heat_capacity > 0) {
                heat_capacity_self_to_sharer += gas_heat_capacity;
            } else {
                heat_capacity_sharer_to_self -= gas_heat_capacity;
            }
        }
    }
    if(!immutable) moles = moles_copy;
    if(!sharer.immutable) sharer.moles = sharer_copy;
	last_share = abs_moved_moles;

    if(abs_temperature_delta > MINIMUM_TEMPERATURE_DELTA_TO_CONSIDER) {
        float new_self_heat_capacity = old_self_heat_capacity + heat_capacity_sharer_to_self - heat_capacity_self_to_sharer;
		float new_sharer_heat_capacity = old_sharer_heat_capacity + heat_capacity_self_to_sharer - heat_capacity_sharer_to_self;

		//transfer of thermal energy (via changed heat capacity) between self and sharer
		if(!immutable && new_self_heat_capacity > MINIMUM_HEAT_CAPACITY) {
			temperature = (old_self_heat_capacity*temperature - heat_capacity_self_to_sharer*temperature_archived + heat_capacity_sharer_to_self*sharer.temperature_archived)/new_self_heat_capacity;
		}

		if(!sharer.immutable && new_sharer_heat_capacity > MINIMUM_HEAT_CAPACITY) {
			sharer.temperature = (old_sharer_heat_capacity*sharer.temperature-heat_capacity_sharer_to_self*sharer.temperature_archived + heat_capacity_self_to_sharer*temperature_archived)/new_sharer_heat_capacity;
		}
		//thermal energy of the system (self and sharer) is unchanged

		if(std::abs(old_sharer_heat_capacity) > MINIMUM_HEAT_CAPACITY) {
			if(std::abs(new_sharer_heat_capacity/old_sharer_heat_capacity - 1) < 0.1) { // <10% change in sharer heat capacity
				temperature_share(sharer, OPEN_HEAT_TRANSFER_COEFFICIENT);
			}
		}
    }
    garbage_collect();
    sharer.garbage_collect();
    if(temperature_delta > MINIMUM_TEMPERATURE_TO_MOVE || std::abs(moved_moles) > MINIMUM_MOLES_DELTA_TO_MOVE) {
		float our_moles = total_moles();
		float their_moles = sharer.total_moles();
		return (temperature_archived*(our_moles + moved_moles) - sharer.temperature_archived*(their_moles - moved_moles)) * R_IDEAL_GAS_EQUATION / volume;
    }
    return 0;
}

float GasMixture::temperature_share(GasMixture &sharer, float conduction_coefficient) {
    float temperature_delta = temperature_archived - sharer.temperature_archived;
    if(std::abs(temperature_delta) > MINIMUM_TEMPERATURE_DELTA_TO_CONSIDER) {
        float self_heat_capacity = heat_capacity_archived();
        float sharer_heat_capacity = sharer.heat_capacity_archived();

        if((sharer_heat_capacity > MINIMUM_HEAT_CAPACITY) && (self_heat_capacity > MINIMUM_HEAT_CAPACITY)) {
            float heat = conduction_coefficient * temperature_delta * (self_heat_capacity*sharer_heat_capacity/(self_heat_capacity+sharer_heat_capacity));
            if(!immutable)
                temperature = std::max(temperature - heat/self_heat_capacity, TCMB);
            if(!sharer.immutable)
                sharer.temperature = std::max(sharer.temperature + heat/sharer_heat_capacity, TCMB);
        }
    }
    return sharer.temperature;
}

float GasMixture::temperature_share(float conduction_coefficient,float sharer_temperature,float sharer_heat_capacity)
{
    float temperature_delta = temperature_archived - sharer_temperature;
    if(std::abs(temperature_delta) > MINIMUM_TEMPERATURE_DELTA_TO_CONSIDER) {
        float self_heat_capacity = heat_capacity_archived();

        if((sharer_heat_capacity > MINIMUM_HEAT_CAPACITY) && (self_heat_capacity > MINIMUM_HEAT_CAPACITY)) {
            float heat = conduction_coefficient * temperature_delta * (self_heat_capacity*sharer_heat_capacity/(self_heat_capacity+sharer_heat_capacity));
            if(!immutable)
                temperature = std::max(temperature - heat/self_heat_capacity, TCMB);
            sharer_temperature = std::max(sharer_temperature + heat/sharer_heat_capacity, TCMB);
        }
    }
    return sharer_temperature;
}

int GasMixture::compare(GasMixture &sample) const {
	float our_moles = 0;
	for (int i = 0; i < moles.size(); i++) {
		float gas_moles = moles[i];
		float delta = std::abs(gas_moles - sample.moles[i]);
		if (delta > MINIMUM_MOLES_DELTA_TO_MOVE && (delta > gas_moles * MINIMUM_AIR_RATIO_TO_MOVE)) {
			return i;
		}
		our_moles += gas_moles;
	}
	if (our_moles > MINIMUM_MOLES_DELTA_TO_MOVE) {
		float temp_delta = std::abs(temperature - sample.temperature);
		if (temp_delta > MINIMUM_TEMPERATURE_DELTA_TO_SUSPEND) {
			return -1;
		}
	}
	return -2;
}

void GasMixture::clear() {
	if (immutable) return;
	moles.clear();
    moles_archived.clear();
}

void GasMixture::multiply(float multiplier) {
	if (immutable) return;
    std::transform(std::execution::seq,
        moles.begin(),moles.end(),
        moles.begin(),
        [&multiplier](auto& gas) {
            return gas*multiplier;
        });
}