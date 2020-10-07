#pragma once

#include "GasMixture.h"

#include "../core/core.h"

#include "../core/proc_management.h"

#include "../third_party/robin_hood.h"

extern std::vector<Value> gas_id_to_type;

extern std::unordered_map<unsigned int, int> gas_ids;

extern int total_num_gases;

class Reaction
{
    public:
        bool check_conditions(const GasMixture& mix) const;
        int react(GasMixture& mix,Value src,Value holder) const;
        inline float get_priority() { return priority; }
        std::string get_name() { return name; };
        Reaction(Value v);
    private:
        Reaction();
        std::string name;
        unsigned int major_gas;
        float priority;
        float min_temp_req = 0.0;
        float max_temp_req = std::numeric_limits<float>::max();
        float min_ener_req = 0.0;
        robin_hood::unordered_map<int,float> min_gas_reqs;
        unsigned int proc_id;
};

extern std::vector<Reaction> cached_reactions;