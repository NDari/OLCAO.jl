import FortranFiles: FortranFile, read

#= struct SpaceGroup =#
#=     space_lattice::Char =#
#=     space_group_number::Integer =#
#=     num_space_ops::Integer =#
#=     num_shifts::Integer =#
#=     space_ops::Array{Float64, 3} =#
#=     space_shifts::Array{Float64, 2} =#
#= end =#

function main(filename::String)
    f = FortranFile(filename)

    # subroutine readSpaceOps
    space_lattice = read(f, Char)
    space_group_number = read(f, Integer)
    num_space_ops, num_shifts = read(f, Integer, Integer)

    space_ops = Array{Float64, 3}(undef, 3, 3, num_space_ops)
    space_shifts = Array{Float64, 2}(undef, 3, num_space_ops)

    for i = 1:num_space_ops
        read(f)
        space_ops[:, 1, i] = read(f, (Float64, 3))
        space_ops[:, 2, i] = read(f, (Float64, 3))
        space_ops[:, 3, i] = read(f, (Float64, 3))
        space_shifts[:, i] = read(f, (Float64, 3))
    end

    # subroutine readLattice
    make_full = read(f, Integer)
    real_lattice = read(f, (Float64, 3, 3))

    # subroutine readAtomicData
    num_atoms = read(f, Integer)

    atom_fract_abc = Array{Float64, 2}(undef, 3, num_atoms)
    atom_element_id = Array{Float64, 1}(undef, num_atoms)
    atom_species_id = Array{Float64, 1}(undef, num_atoms)

    for i = 1:num_atoms
        atom_element_id[i], atom_species_id[i], atom_fract_abc[:, i] = read(f, (Float64, 5))
    end

    # subroutine prepareSymmAtomicCoords
    symm_fract_abc = Array{Float64, 2}(undef, 3, num_atoms*num_space_ops)
    symm_element_id = Array{Float64, 1}(undef, num_atoms*num_space_ops)
    symm_species_id = Array{Float64, 1}(undef, num_atoms*num_space_ops)

    num_symm_atoms = 0

    # subroutine applySymmetry
    new_position = zeros(Float64, 3)
    for i = 1:num_atoms
        for j = 1:num_space_ops
            for k = 1:3
                new_position[k] = zero(Float64)

                for l = 1:3
                    new_position[k] += atom_fract_abc[l, i] * space_ops[l, k, j]
                end

                new_position[k] += space_shifts[k, j]
            end

            shift_to_cell(new_position)

            if find_match(new_position, num_symm_atoms, symm_fract_abc) == 0
                num_symm_atoms += 1
                symm_fract_abc[:, num_symm_atoms] = new_position[:]
                symm_element_id[num_symm_atoms] = atom_element_id[i]
                symm_species_id[num_symm_atoms] = atom_species_id[i]
            end
        end

        # subroutine reduceCell

    end



    #= return SpaceGroup(space_lattice, =#
    #=                   space_group_number, =#
    #=                   num_space_ops, =#
    #=                   num_shifts, =#
    #=                   space_ops, =#
    #=                   space_shifts) =#
end

function shift_to_cell(new_position, small_thrash = 1.0e-8)
    for i = 1:3
        if abs(new_position[i] < small_thrash)
            new_position[i] = 0.0
        elseif abs(new_position[i] - 1.0) < small_thrash
            new_position[i] = 0.0
        end
    end

    for i = 1:3
        if new_position[i] < 0.0
            new_position[i] += 1.0
        elseif  new_position[i] > 1.0
            new_position[i] -= 1.0
        end
    end
end

function find_match(new_position, num_symm_atoms, symm_fract_abc, small_thresh = 1.0e-8)
    for i = 1:num_symm_atoms
        if sum((new_position[:] - symm_fract_abc[:, i])^2) < small_thresh
            return 1
        end
    end
    return 0
end

function reduce_cell(make_full, space_lattice)
    if make_full == 1 || space_lattice == 'P' || space_lattice == 'H'
        return
    end
end
