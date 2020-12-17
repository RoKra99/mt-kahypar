#pragma once

#include "mt-kahypar/definitions.h"


namespace mt_kahypar {
namespace ds {

template <typename HashMap>
class CommunityCount {

public:
    CommunityCount();
    ~CommunityCount();


private:

    std::vector<PartitionID> _single_cut_communities;
    HashMap _communities;

};

}
}