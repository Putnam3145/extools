#pragma once

#include "GasMixture.h"

#include "../core/core.h"

#include "../core/proc_management.h"

extern std::vector<Value> gas_id_to_type;

extern int total_num_gases;

class Reaction
{
    public:
        bool check_conditions(const GasMixture& mix) const;
        int react(GasMixture& mix,Value src,Value holder) const;
        inline float get_priority() { return priority; }
        Reaction(Value v);
    private:
        Reaction();
        int major_gas;
        float priority;
        float min_temp_req = 0.0;
        float max_temp_req = std::numeric_limits<float>::max();
        float min_ener_req = 0.0;
        robin_hood::unordered_map<int,float> min_gas_reqs;
        unsigned int proc_id;
};

extern std::vector<Reaction> cached_reactions;