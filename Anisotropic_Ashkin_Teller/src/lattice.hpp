/** 
 * @brief class constructing a regular lattice with arbitrary dimensions
 *
 * @author Marco Mueller
 * @date 11.04.2011
 *
 */

# ifndef __LATTICE_HPP_
# define __LATTICE_HPP_

# include "const.hpp"
# include "simulation.hpp"

class lattice
{
    private:

        // all information about physics, geometry, ...
        // are stored in "parameters"
        simulation parameters;

        // the spins of the lattice are mapped onto an one-dimensional 
        // vector of integers called "state"
//        int32_1d state;
        // to abstract the lattice geometry neighbour tables are used
        // that store the indices of the neighbours in the state - vector

        // vec[spin] = vector of (diagonal) neighbours
        uint32_2d near_neighbours; ///< nearest neighbours

        // fixed boundary conditions are implemented by widening the hypercube
        // in each dimension:
        uint32 correct_boundary_index(int32_1d &, uint32);

        uint32_1d dimensions;
        uint32 dimension;
        uint32 volume;

        template<class T> int32 cartesian2index(T &ccord);
        uint32 get_virtual_index(int32_1d &ccord);

    public:
        // this constructor calculates the lattice geometry
        lattice(simulation &parameters);
        void calculate_neighbours();
        void print_grid_list(std::ostream &out);
        
        inline uint32 get_near_neighbour(uint32 site, uint32 neighbour);

        inline uint32 number_of_near_neighbours();
        
        uint32 number_nn;

};

/// Constructor.
/// Constructs the neighbour tables according to the boundary conditions
/// and dimensions.
lattice::lattice(simulation &parameters_) : parameters(parameters_)
{
    // ask for lattice dimensions
    parameters.get_dimensions(dimensions);

    // calculate the dimension of the lattice
    dimension = static_cast<uint32>(dimensions.size());
    volume = parameters.get_number_of_spins();

    calculate_neighbours();
    
    number_nn = static_cast<uint32>(near_neighbours.at(0).size());
}

void lattice::calculate_neighbours()
{

    std::string dimensions_string(" ");
    for (size_t i = 0; i < dimensions.size(); i++)
    {
        dimensions_string += dlib::cast_to_string(dimensions.at(i));
        if (i < dimensions.size()-1) dimensions_string += ", ";
    }

    // initialize neighbour tables
    // convention is vector[sites]([plaquettes])[neighbours]
    //
    // there are 2*d nearest neighbours for each spin
    near_neighbours.resize(volume);

    // Constructs a lattice with fixed boundary conditions 
    // that forces a system to develop interfaces.
    // In the gonihedric ising model this means that the upper half
    // of the hypercube is fixed with +1-spins      + + + + + + + +
    // and the lower half has -1-spins.             + o o o o o o +
    // This also means that at least in one         + o o o o o o +
    // direction the number of spins has to be      - o o o o o o -
    // even. The graphic on the right hand side     - o o o o o o -
    // shows the interface boundary condition for   - - - - - - - -
    // a two dimensional grid with dimensions 4x6.

    // Constructs a lattice with periodic boundary conditions
    
    // stores cartesian coordinates of the lattice spin
    int32_1d ccord(dimension,0);

    for (uint32 site=0; site < volume; site++)
    {
        // volume of the current slice under consideration
        int32 dim_vol = 1;
        // iterate over all dimensions
        for (uint32 q = 0; q < dimension; ++q)
        {
            // check whether the boundary was hit
            if (site % dim_vol == 0 && site > 0)
            {
                ccord[q]++;
                if (q > 0)
                    ccord[q-1] = 0;
            }
            dim_vol *= dimensions.at(q);
        }
        
        
        for (uint32 q = 0; q < dimension; ++q)
        {
            // cartesian coordinates of the left and right neighbour
            // initialized by the coordinate of 
            // the site under consideration
            int32_1d ccord_left;
                ccord_left.assign(ccord.begin(),ccord.end());
            int32_1d ccord_right;
                ccord_right.assign(ccord.begin(),ccord.end());

            // ---
            // nearest neighbours:
            // ---
            // increase (right neighbours) or decrease (left neighbours)
 
            ccord_left[q] = 
                ccord_left[q] - 1;
            ccord_right[q] = 
                ccord_right[q] + 1;

            // now calculate the index of the neighbours
            uint32 left = cartesian2index(ccord_left);
            uint32 right = cartesian2index(ccord_right);
            
            // check for boundary conditions 
            if (ccord_left[q] < 0)
                left  = correct_boundary_index(ccord_left, q);
            if (ccord_right[q] >= static_cast<int32>(dimensions.at(q)))
                right = correct_boundary_index(ccord_right, q);

            // store the calculated indices
            near_neighbours[site].push_back(right);
            near_neighbours[site].push_back(left);
            
        } // for (dimensions)
    } // for (sites)
}


// calculates the index in a 1d - Array when cartesian coordinates are given.
template<class T> int32 lattice::cartesian2index(T &ccord)
{
    int32 index = 0;
    int32 last_dim = 1;
    for (uint32 d = 0; d < dimension; ++d)
    {
        index += last_dim * ccord.at(d);
        last_dim *= static_cast<int32>(dimensions.at(d));
    }
    return index; 
}

uint32 lattice::correct_boundary_index(int32_1d &cc, uint32 q)
{
    if ((parameters.get_boundary_conditions() == INTERFACE)
       || (parameters.get_boundary_conditions() == VANISH))
    {
        return get_virtual_index(cc);
    }
    else if (parameters.get_boundary_conditions() == PERIODIC)
    {
        // Constructs a lattice with periodic boundary conditions
        if (cc[q] < 0) // method was called for a spin on the left boundary
            cc[q] = ( cc.at(q) + dimensions.at(q) ) % dimensions.at(q);
        else // method was called for a spin on a right boundary
            cc[q] = cc.at(q) % dimensions.at(q);
    }
    else if (parameters.get_boundary_conditions() == HELICAL)
    {
        const int ScrewParameter = 1;
        int dimScrew = (q - 1 + dimensions.size()) % dimensions.size();
        if (cc[q] < 0) // method was called for a spin on the left boundary
        {
            cc[q] = ( cc.at(q) + dimensions.at(q) ) % dimensions.at(q);
            cc.at(dimScrew) = (cc.at(dimScrew) + ScrewParameter) 
                % dimensions.at(dimScrew);
        }
        else // method was called for a spin on the right boundary
        {
            cc[q] =  cc.at(q) % dimensions.at(q);
            cc.at(dimScrew) = 
                (cc.at(dimScrew) - ScrewParameter + dimensions.at(dimScrew)) 
                    % dimensions.at(dimScrew);
        }
    }
    return cartesian2index(cc);
}

// returns the virtual index of a site at cartesian coordinates ccord
uint32 lattice::get_virtual_index(int32_1d &ccord)
{
    if (ccord[0] < static_cast<int32>(dimensions[0]/2))
        return volume;
    else 
        return volume + 1;
}

///prints a table of spinstates and neighbours
void lattice::print_grid_list(std::ostream &out)
{
	out << "# index, near n.; next near n.\n";
    // uint32 dimension = parameters.get_dimension();
	for (uint32 i=0; i < parameters.get_number_of_spins(); ++i)
	{
		out << i << ", ";
		for (uint32 q = 0; q < dimension; ++q)
        {
			out << near_neighbours.at(i).at(q) << ", ";
            out << near_neighbours.at(i).at(dimension+q);
            (q < dimension -1)? out << ", " : out << "; ";
		}
        out << "\n";
	}
}

inline uint32 lattice::get_near_neighbour(uint32 site, uint32 neighbour)
{
    return near_neighbours[site][neighbour];
}

inline uint32 lattice::number_of_near_neighbours()
{
    return number_nn;
}

# endif // __LATTICE_HPP__
