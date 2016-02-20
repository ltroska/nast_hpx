#ifndef UTIL_HELPERS_HPP
#define UTIL_HELPERS_HPP

#define FOR_EVERY_BOUNDARY_PARTITION    for (uint l = 1; l < p.num_partitions_y - 1; l++) \
        for (uint k = 1; k < p.num_partitions_x - 1; k++) \
        { \
            if (!(k == 1 || k == p.num_partitions_x - 2 || l == 1 || l == p.num_partitions_y -1))\
                continue;

#define FOR_EVERY_PARTITION     for (uint l = 1; l < p.num_partitions_y - 1; l++) \
        for (uint k = 1; k < p.num_partitions_x - 1; k++)

#define END_FOR }

//needed for deleting the partition_data arrays
template<typename T>
class array_deleter {
    public:
        void operator()(T const* p)
        {
            delete [] p;
        }
};


//check if point is in global id range
inline bool in_range(uint start_i, uint end_i, uint start_j, uint end_j, uint i,  uint j)
{
    if (i < start_i || i > end_i || j < start_j || j > end_j)
        return false;

    return true;
}

//finding id of neighboring localities
inline uint get_neighbor_id(uint id, direction dir, uint num_localities)
{
    uint res_x, res_y;
    if (num_localities == 2)
    {
        res_x = 2;
        res_y = 1;
    }
    else
    {
        res_x = static_cast<uint>(sqrt(num_localities));
        res_y = res_x;
    }

    switch (dir)
    {
        case LEFT:
            return ((id-1)/res_x == id/res_x && id-1 < num_localities) ? id-1 : num_localities;

        case RIGHT:
            return ((id+1)/res_x == id/res_x && id+1 < num_localities) ? id+1 : num_localities;

        case TOP:
            return ((id+res_x) < num_localities && (id+res_x)/res_x == id/res_x +1) ? id+res_x : num_localities;

        case BOTTOM:
            return ((id-res_x) < num_localities && (id-res_x)/res_x == id/res_x -1) ? id-res_x : num_localities;

        case TOP_LEFT:
            return ((id+res_x-1) < num_localities && (id+res_x-1)/res_x == id/res_x +1) ? id+res_x-1 : num_localities;

        case TOP_RIGHT:
            return ((id+res_x+1) < num_localities && (id+res_x+1)/res_x == id/res_x+1) ? id+res_x+1 : num_localities;

        case BOTTOM_LEFT:
            return ((id-res_x-1) < num_localities && (id-res_x-1)/res_x == id/res_x-1) ? id-res_x-1 : num_localities;

        case BOTTOM_RIGHT:
            return ((id-res_x+1) < num_localities && (id-res_x+1)/res_x == id/res_x-1) ? id-res_x+1 : num_localities;
        default:
            return num_localities;
    }
}
#endif
