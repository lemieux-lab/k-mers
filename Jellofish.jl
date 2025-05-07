using Pkg
Pkg.activate("k_mers/k_mers/")

using JuBox
using BioSequences
using BioSymbols
using Kmers
using ProgressMeter
using Base.Threads

struct Node{Sb<:BioSymbol}
    value::Sb
    nexts::Vector{Node{Sb}}
    count_idx::UInt32
end

struct Trie{Sb<:BioSymbol}
    K::Int
    roots::Vector{Node{Sb}}
    counts::Vector{UInt32}
end

function Base.push!(trie::Trie{Sb}, kmer::Km) where {Sb<:BioSymbol, Km<:Kmer}
    nexts = trie.roots
    last = nothing
    for (i, n) in enumerate(kmer)
        next = get_next(nexts, n)
        # println(next)
        if isnothing(next)  # Node does not yet exist
            if i == length(kmer)  # Node is a leaf, need to add a count entry
                push!(trie.counts, 0)
                # println(trie.counts)
                next = Node{Sb}(n, Node{Sb}[], length(trie.counts))
            else
                next = Node{Sb}(n, Node{Sb}[], 0)
            end
            push!(nexts, next)
        end
        last = next
        nexts = next.nexts
    end
    trie.counts[last.count_idx] += 1
    
end

function get_next(v::Vector{Node{Sb}}, symbol::Sb) where {Sb<:BioSymbol}
    for n in v
        n.value == symbol && return n
    end
end

function init_trie(K::Int, Sb::Type{<:BioSymbol}=DNA)
    alpha = Node{Sb}[Node{Sb}(symbol, Node{Sb}[], 0) for symbol in alphabet(Sb) if !(isambiguous(symbol) || isgap(symbol))]
    return Trie{Sb}(K, alpha, UInt32[])
end

function k_merize(sequence::LongSequence{Ab}; K::Int) where {Ab<:Alphabet}
    length(sequence) < K && error("Read of size $(length(sequence)) too small for $K-mers")
    to_return = Kmer{Ab, K}[]
    for i in 1:(length(sequence)-K+1)
        push!(to_return, Kmer{Ab, K}(sequence[i:i+K-1]))
    end
    return to_return
end

function jello_trie(fastq::String, Sb::Type{<:BioSymbol}=DNA; K::Int)
    GC.enable(false)
    trie = init_trie(K, Sb)
    wcl = countlines(fastq)
    progess = Progress(wcl, desc="Parsing fastq into trie...")
    open(fastq) do f
        for (i, l) in enumerate(eachline(f))
            (i-2)%4 != 0 && continue    # Skips non-sequence lines of fastq
            
            seq = bioseq(l)
            kmers = k_merize(seq, K=K)
            push!.(Ref(trie), kmers)
            next!(progess; showvalues=[
                ("unique k-mers", length(trie.counts))
                ("counts memory size", Base.format_bytes(Base.summarysize(trie.counts)))
                ])
        end
    end
    GC.enable(true)
    return trie
end

function collapse(kct::Kct)

end

function concurrent_merger(merge_queue::Channel{Kct{1, K}}) where K
    queue = Vector{Kct{1, K}}()
    lock = ReentrantLock()
    cond = Condition()

    # Consumer that listens to the input channel and accumulates KCTs
    reader = @async begin
        for kct in merge_queue
            lock(lock) do
                push!(queue, kct)
                notify(cond)
            end
        end
    end

    # Parallel merger tasks
    merger = @async begin
        while true
            local_kct1, local_kct2 = nothing, nothing

            lock(lock) do
                while length(queue) < 2 && isopen(merge_queue)
                    wait(cond)
                end
                if length(queue) >= 2
                    local_kct1 = pop!(queue)
                    local_kct2 = pop!(queue)
                elseif !isopen(merge_queue) && length(queue) == 1
                    return queue[1]  # final mother_kct
                end
            end

            if local_kct1 !== nothing && local_kct2 !== nothing
                @async begin
                    merged = merge(local_kct1, local_kct2)
                    lock(lock) do
                        push!(queue, merged)
                        notify(cond)
                    end
                end
            end
        end
    end

    return fetch(merger)  # final mother_kct
end


function jello_kct_chunk(fastq::String; K::Int, chunking::Int=1_000_000, queue_size::Int=32)
    lines = readlines(open(fastq, "r"))[2:4:end]
    chunks = [lines[i:min(i+chunking-1, end)] for i in 1:chunking:length(lines)]
    merge_queue = Channel{Kct{1, K}}(queue_size)
    mother_kct = Kct(1, K)
    # Merger thread. Works asynchronously from workers.
    merger = @async try
        for child_kct in merge_queue
            # error("TEST TEST TEST")
            # println("merging child with k_mers: $(length(child_kct.table))")
            mother_kct = merge(mother_kct, child_kct)
        end
    catch e
        @error "Merger task crashed" exception=(e, catch_backtrace())
    end

    progress = Progress(length(chunks), desc = "Processing $chunking k-mers chunks...")
    # Worker threads, building child Kcts asynchronously
    @threads for chunk in chunks
        child_kct = Kct(1, K)
        for l in chunk
            # println(l)
            seq = bioseq(l)
            DNA_N in seq && continue
            k_mers = k_merize(seq, K=K)
            for k_mer in k_mers
                buf = JuBox.JaggedRLEArray(UInt32)
                push!(buf, UInt32[1])
                id = push!(child_kct, buf)
                # println(k_mer)
                push!(child_kct.table, JuBox.KRecord(k_mer, id))
            end
        end
        # println("starting a sort")
        sort!(child_kct)
        # println(merge_queue.n_avail_items)
        put!(merge_queue, child_kct)
        # println("finished one chunk")
        next!(progress; showvalues=[
                ("items in merge queue", merge_queue.n_avail_items)
                ])
    end

    close(merge_queue)
    wait(merger)

    return mother_kct
end