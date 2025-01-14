# New version of KCT supporting new Kmers implementation, variable kmer size, variable kmer type, (RLE?)
using Pkg
Pkg.activate("k_mers/k_mers/")

using Kmers
using BioSequences
using JSON
using ElasticArrays
using ProgressMeter
include("parallel_sort.jl")

export Kct, Kct, load, save, merge
const seq_codes = Dict{String, Type}(
    "AA"=>AminoAcidAlphabet,
    "DNA"=>DNAAlphabet{2},
    "DNA4"=>DNAAlphabet{4},
    "RNA"=>RNAAlphabet{2},
    "RNA4"=>RNAAlphabet{4},
)

struct KRecord{K, Ab}
    seq::Mer{K, Ab}
    pool_id::UInt32
end

Base.isless(a::KRecord{K, Ab}, b::KRecord{K, Ab}) where {K, Ab<:Type} = a.seq < b.seq 
Base.isless(a::Mer{K, Ab}, b::KRecord{K, Ab}) where {K, Ab<:Type} = a < b.seq
Base.isless(a::KRecord{K, Ab}, b::Mer{K, Ab}) where {K, Ab<:Type} = a.seq < b
Base.:(==)(a::KRecord{K, Ab}, b::KRecord{K, Ab}) where {K, Ab<:Type} = a.seq == b.seq

# const mask_big = 0x00000001 << 27
isbig(x::UInt32, mask_big::UInt32) = x & mask_big != 0
makebig(x::Integer, mask_big::UInt32) = UInt32(x) | mask_big
idbig(x::UInt32, mask_big::UInt32) = x & (~mask_big)

struct Kct{N, K, Ab<:Type}
    table::Vector{KRecord{K, Ab}}
    idx::Vector{UnitRange{Int64}}
    small::ElasticMatrix{UInt8}
    big::ElasticMatrix{UInt32}
end

Kct(N, K, A::String) = Kct(N, K, seq_codes[A])
Kct(N, K, A::Ab=DNAAlphabet{2}, t=0, s=0, b=0) where {Ab<:Type} = Kct{N, K, Ab}(
    Vector{KRecord{K, A}}(undef, t),
    fill(0:-1, 4^11),
    Matrix{UInt8}(undef, N, s),
    Matrix{UInt32}(undef, N, b)
)

function Base.push!(kct::Kct{N, K, Ab}, count::Matrix{UInt32}) where {N, K, Ab<:Type}
    if max(count...) ≤ typemax(UInt8)
        append!(kct.small, count)
        (_, id) = UInt32.(size(kct.small))
    else
        append!(kct.big, count)
        (_, id) = UInt32.(size(kct.big))
        id |=  0x00000001 << K
    end
    return id
end

function Base.push!(kct::Kct{N, K, Ab}, count::Matrix{UInt8}) where {N, K, Ab<:Type}
    append!(kct.small, count)
    (_, id) = UInt32.(size(kct.small))
    return id
end

function Base.push!(kct::Kct{N, K, Ab}, seq::Mer{K, Ab}, count::Matrix{UInt32}) where {N, K, Ab<:Type}
    id = push!(kct, count)
    push!(kct.table, KRecord(seq, id))
    return id
end

function Base.getindex(kct::Kct{N, K, Ab}, i::Integer) where {N, K, Ab<:Type}
    rec = kct.table[i]
    if mask_big & rec.pool_id != 0
        tmp_t = convert(Vector{UInt32}, kct.big[:, rec.pool_id & (~mask_big)])
    else
        tmp_t = convert(Vector{UInt32}, kct.small[:, rec.pool_id])
    end
    return rec.seq, tmp_t
end


function Base.getindex(kct::Kct{N, K, Ab}, idx::AbstractArray{M}) where {N, K, Ab<:Type, M<:Integer}
    return [kct[i] for i in idx]
end

#= 
Deletes an element in the table, but does not check if the corresponding count vector
is no longer being pointed at by any element in the table. This can lead to garbage 
data left in memory until the whole Kct is freed.
Do not rely on this to make a Kct smaller. This is just a convenience thing
for when we need to work with indexes on a Kct.
=#
function Base.deleteat!(kct::Kct{N, K}, idx::Int64) where {N, K}
    deleteat!(kct.table, idx)
end

function Base.deleteat!(kct::Kct{N, K}, idx::AbstractArray{M}) where {N, K, M<:Integer}
    for i in idx
        deleteat!(kct, i)
    end
end

#= 
deleteat scales very powerly with big kcts, since it requires all elements beyond the
deleted one to be moved in memory.
quick_delete instead moves the last element to the position of the elmenent to delete
and then cuts of the tail of the table, meaning no re-assignation nescessary, which
means this is in constant time, regardless of kct size (which is pretty nice).
However, this does NOT preserve order in the table, which means this should
never be used before the kct is finished merging with others, or if you
want to be able to search k-mers using dichotomic search.
=#
function quick_delete!(kct::Kct{N, K}, idx::Int64) where {N, K}
    kct.table[idx] = kct.table[end]
    deleteat!(kct.table, length(kct.table))
end

function quick_delete!(kct::Kct{N, K}, idx::AbstractArray{M}) where {N, K, M<:Integer}
    for i in idx
        quick_delete!(kct, i)
    end
end

function Base.iterate(kct::Kct, state=1)
    state > length(kct) && return nothing
    return (kct[state], state+1)
end

Base.length(kct::Kct) = length(kct.table)

# Base.getindex(kct::Kct{N,T}, v::Vector{U}) where {N,T,U<:Integer} = [kct[i] for i in v]

# function Base.getindex(kct::Kct{N, T}, key::DNAKmer{K, 1}{K})::Tuple{DNAKmer{K, 1}{K}, NTuple{N, UInt32}} where {N, T}
#     idx_key = (key.data[1] >> 40) + 1
#     r = kct.idx[idx_key]
#     t = @view kct.table[r]
#     i = searchsortedfirst(t, key)
#     if i > length(t) || t[i].seq != key
#         return (DNAKmer{K, 1}{K}("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), NTuple{N, UInt32}(zeros(UInt32, N)))
#     end
#     first(kct.pool[t[i].pool_id]) == typemax(T) && return kct.next[key]
#     return (t[i].seq, convert(NTuple{N, UInt32}, kct.pool[t[i].pool_id]))
# end

# Base.getindex(kct::Kct, key::String) = kct[DNAKmer{K, 1}{K}(key)]

function Base.show(io::IO, t::Tuple{Mer{K, Ab}, NTuple{N, UInt32}}) where {N, K, Ab<:Type}
    print(io, t[1], ": ", join(string.(t[2]), ", "))

end

function Base.sort!(kct::Kct)
    psort!(kct.table)
end

function computeIndex!(kct::Kct{N, K, Ab}) where {N, K, Ab<:Type}
    start = 1
    last_key = 0x0000000000000000
    bits_per_symbol = BioSequences.bits_per_symbol(BioSequences.BitsPerSymbol(Ab()))
    
    for i=1:length(kct.table)
        key = kct.table[i].seq.data[1] >> 20*bits_per_symbol
        if key > last_key
            kct.idx[last_key+1] = start:i
            start = i
            last_key = key
        end
    end
end

function _dedup(kct, nb)
    o = psortperm(nb)
    last = 0
    id = 0
    buf = Matrix{UInt32}(undef, 1, 1)
    for i in o
        if nb[i] != last
            buf[1, 1] = nb[i]
            id = push!(kct, buf)
        end
        kct.table[i] = KRecord(kct.table[i].seq, id)
        last = nb[i]
    end

    return kct
end

"""
    Kct(fn::String)

Parse a Jellyfish K-mer count table into a sorted array of DNAKmer{K, 1}. No
verification is made on the header, use with caution. 
"""
function Kct(fn::String, A::Ab=DNAAlphabet{2}; big_only::Bool=false) where {Ab<:Type}
    
    f = open(fn, "r")
    offset = parse(Int, readuntil(f, "{"))
    json_start = position(f)
    seek(f, position(f) - 1) # Need to grab the {
    header = JSON.parse(f)
    println(header)
    seek(f, offset + json_start - 1) 

    K = parse(Int, header["cmdline"][findfirst(x->x=="-m", header["cmdline"])+1])
    out_counter_len = findfirst(x->x=="--out-counter-len", header["cmdline"])
    count_bytes = isnothing(out_counter_len) ? 4 : parse(Int, header["cmdline"][out_counter_len+1])
    println("Found k-mers of size $K and count encoded on $count_bytes bytes.")

    count_type = Dict(1=>UInt8,
                  2=>UInt16,
                  4=>UInt32)[count_bytes]

    count_mask = Dict(1=>0xFF,
                  2=>0xFFFF,
                  4=>0xFFFFFFFF)[count_bytes]

    count_shift = 64 - count_bytes * 8

    nb = Vector{count_type}()
    kct = Kct(1, K)
    prog = ProgressUnknown()

    # Kmer{DNAAlphabet{2}, 31, 1}(Kmers.unsafe, NTuple{1, UInt64}(0x0000000000000000,))
    # N term of kmers is the number of values in the Tuple (DNA: K/32, AA: K/8 for Int64)
    while(!eof(f))
        if count_shift - K*2 < 0
            println((read(f, UInt64),))
            tmp_seq = Kmer{A, K, N}((read(f, UInt64),))
            tmp_nb = read(f, count_type)
        else
            line_bits = (read(f, UInt64),)
            tmp_seq = Mer{K, A}(line_bits)
            tmp_nb = count_type(line_bits[1] >> count_shift & count_mask)
        end

        push!(kct.table, KRecord(tmp_seq, UInt32(0)))
        push!(nb, tmp_nb)
        next!(prog)
    end

    _dedup(kct, nb)
    
    sort!(kct)

    # computeIndex!(kct)

    return kct
end

function Base.write(s::IO, kct::Kct{N, K}) where {N, K}
    write(s, N)
    write(s, length(kct.table))
    write(s, size(kct.small)[2])
    write(s, size(kct.big)[2])
    write(s, kct.table)
    write(s, kct.small)
    write(s, kct.big)
end

function Base.read(s::IO, ::Type{Kct})
    N = read(s, Int)
    l = read(s, Int)
    ls = read(s, Int)
    lb = read(s, Int)
    kct = Kct(N, l, ls, lb)
    read!(s, kct.table)
    read!(s, kct.small)
    read!(s, kct.big)
    computeIndex!(kct)
    return kct
end

function Base.findfirst(kct::Kct{N, K}, key::Mer{K, Ab}) where {N, K, Ab<:Type}
    idx_key = (key.data[1] >> 40) + 1
    r = kct.idx[idx_key]
    t = @view kct.table[r]
    i = searchsortedfirst(t, key)
    if i > length(t) || t[i].seq != key
        return 0
    end
    return r[1] + i-1
end

"""
    save(kct::Kct, fn::String)

Saves a K-mer count table into a binary format file.
"""
function save(kct::Kct{N}, fn::String) where {N}
    f = open(fn, "w")
    write(f, kct)
    close(f)
end

"""
    load(fn::String)::Kct

Loads a K-mer count table from a binary format file.
"""
function load(fn::String)::Kct
    f = open(fn, "r")
    return read(f, Kct)
end

# function _getpool(kct::Kct, pool_id::UInt32)::Union{ElasticVector{UInt8}, ElasticVector{UInt32}}
#     isbig(pool_id) && return kct.big[:,idbig(pool_id)]
#     return kct.small[:,pool_id]
# end

function _mergededup(kct::Kct, a::Kct{M, K}, b::Kct{N, K}, p) where {M, N, K}
    o = psortperm(p)
    last::Tuple{UInt32, UInt32} = (0, 0)
    id::UInt32 = 0
    mask_big = 0x00000001 << K
    for i in o
        if p[i] != last
            ia, ib = p[i]
            if isbig(ia, mask_big) || isbig(ib, mask_big)
                bbuf = zeros(UInt32, M+N, 1)
                if isbig(ia, mask_big)
                    bbuf[1:M,1] = a.big[:,idbig(ia, mask_big)]
                elseif ia != 0
                    bbuf[1:M,1] = a.small[:,ia]
                end
                if isbig(ib, mask_big)
                    bbuf[(M+1):(M+N),1] = b.big[:,idbig(ib, mask_big)]
                elseif ib != 0
                    bbuf[(M+1):(M+N),1] = b.small[:,ib]
                end
                id = push!(kct, bbuf)
            else
                buf = zeros(UInt8, M+N, 1)
                if ia != 0
                    buf[1:M,1] = a.small[:,ia]
                end
                if ib != 0
                    buf[(M+1):(M+N),1] = b.small[:,ib]
                end
                id = push!(kct, buf)
            end
        end
        kct.table[i] = KRecord(kct.table[i].seq, id)
        last = p[i]
    end
end

function Base.merge(a::Kct{M, K}, b::Kct{N, K}) where {M, N, K}
    i = 1
    j = 1

    kct = Kct(M+N)
    p = Vector{Tuple{UInt32, UInt32}}()

    while i ≤ length(a) || j ≤ length(b)
        if j > length(b) || (i ≤ length(a) && a.table[i].seq < b.table[j].seq)
            push!(p, (a.table[i].pool_id, 0))
            # K = length(a.table[i].seq)
            push!(kct.table, KRecord(a.table[i].seq, 0))
            i += 1
        elseif i > length(a) || (j ≤ length(b) && b.table[j].seq < a.table[i].seq)
            push!(p, (0, b.table[j].pool_id))
            # K = length(a.table[i].seq)
            push!(kct.table, KRecord(b.table[j].seq, 0))
            j += 1
        else
            push!(p, (a.table[i].pool_id, b.table[j].pool_id))
            # K = a.table[i]
            push!(kct.table, KRecord(a.table[i].seq, 0))
            i += 1
            j += 1
        end
    end

    _mergededup(kct, a, b, p)
    # return kct, p
    return kct
end