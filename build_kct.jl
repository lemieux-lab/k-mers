using Pkg
Pkg.activate("k_mers/k_mers/")

include("Kct.jl")

function build_kct_binary_tree(jf_files::Array{String, 1}, save_path::String="."; big_only::Bool=false)
    samples = [split(path, "/")[9] for path in jf_files]
    done = Set{Int}()

    # Init stuff
    prog = Progress(length(jf_files), "Parsing Jellyfish Files (binary tree version)"); ProgressMeter.update!(prog)

    # Stacks previously loaded Kct
    tree_stack = Array{Kct, 1}()

    # Using a stack to emulate a binary tree merge. Adding a new leaf at every loop
    for (sample_name, file) in zip(samples, jf_files)
        push!(tree_stack, Kct(file, big_only=big_only))
        if !(size(tree_stack[begin].big)[1] in done)
            JuBox.save(tree_stack[begin], save_path*"tcga_brca_$(typeof(tree_stack[begin]))_samples_bt_backup.kct")
            push!(done, size(tree_stack[begin].big)[1])
        end

        # Then checking if that leaf can be merged with another leaf at the same level.
        # Keep merging up until no matching leaf at that level.
        while length(tree_stack) >= 2 && typeof(tree_stack[end]) == typeof(tree_stack[end-1])
            ProgressMeter.update!(prog, showvalues = [("Tree size on RAM: ", Base.format_bytes(Base.summarysize(tree_stack))),
                                        ("Merging", "["*join([typeof(k) for k in tree_stack[1:end-1]], " | ")*" ← $(typeof(tree_stack[end]))]")])
            last = pop!(tree_stack)
            # if !(size(last.big)[1] in done)
            #     JuBox.save(last, save_path*"tcga_brca_$(typeof(last))_samples_bt_backup.kct")
            #     push!(done, size(last.big)[1])
            # end
            push!(tree_stack, (merge(pop!(tree_stack), last)))
        end

        # Updating progress bar
        next!(prog; showvalues = [("Tree size on RAM: ", Base.format_bytes(Base.summarysize(tree_stack)))])
    end
    finish!(prog)
    # Extract the root node
    final_node = popfirst!(tree_stack)
    JuBox.save(final_node, save_path*"tcga_brca_$(typeof(final_node))_samples_bt_backup.kct")
    
    # Cleanup any potential leftover leaf. Merges them linearly regardless of level.
    prog = Progress(length(tree_stack)-1, "Remaining merges to cleanup outstanding leaves"); ProgressMeter.update!(prog)
    for i in 1:length(tree_stack)-1
        pushfirst!(tree_stack, merge(popfirst!(tree_stack), popfirst!(tree_stack)))
        next!(prog)
    end
    finish!(prog)
    JuBox.save(tree_stack[end], save_path*"tcga_brca_$(typeof(tree_stack[end]))_samples_bt_backup.kct")

    # Merging root node and cleanup node, which is the entirely built kct, in order of samples.
    final_node = merge(final_node, pop!(tree_stack))

    # Saving that kct
    JuBox.save(final_node, save_path*"tcga_brca_$(length(samples))_samples_bt.kct")
end

jf_files = readlines("/u/jacquinn/Aboleth/data/multi_tcga_mini_balanced_samples.txt")

build_kct_binary_tree(jf_files, "/u/jacquinn/phd_stuff/data/aboleth_data/multi_tcga_kct_data/new_kct_test_balanced_60")
