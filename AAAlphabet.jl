using Kmers, Kmers.BioSequences

struct AAAlphabet <: Alphabet end

BioSequences.BitsPerSymbol(::AAAlphabet) = BioSequences.BitsPerSymbol{5}()
Base.eltype(::Type{AAAlphabet}) = AminoAcid
bin(x) = Base.bin(x, 64, false)

function BioSequences.encode(::AAAlphabet, aa::AminoAcid)
    reinterpret(UInt8, aa) > reinterpret(UInt8, AA_Gap) && return nothing
    return convert(UInt, reinterpret(UInt8, aa))
end

@inline function BioSequences.decode(::AAAlphabet, x::UInt)
    return reinterpret(AminoAcid, x % UInt8)
end

offset(::Km) where {Km<:Kmer{A, K, N}} where {A<:Alphabet, K, N} = N*64 - BioSequences.bits_per_symbol(A())*K
true_index(k::Km, i::Int) where {Km<:Kmer{A, K, N}} where {A<:Alphabet, K, N} = offset(k)÷BioSequences.bits_per_symbol(A()) + i
carry(c::Int) = c-c÷5
# chunk(k::Km, i::Int) where {Km<:Kmer} = true_index(k, i)÷(12+carry(true_index(k, i)÷12))

function chunk(k::Km, i::Int) where {Km<:Kmer}
    off_idx = true_index(k, i)
    interval = 1:12
    bps = BioSequences.bits_per_symbol(Alphabet(k)) % UInt
    cycle = length(k.data) % bps
    cur_chunk = 0

    while !(off_idx in interval)
        interval = interval[end]+1:interval[end] + 12 + (cycle != 0)
        # println(interval)
        cur_chunk += 1
        cycle = (length(k.data) - cur_chunk) % bps
    end
    # println("\n$interval")
    return cur_chunk
end

@inline function BioSequences.extract_encoded_element(seq::Kmer{<:AAAlphabet}, i::Integer)
    off = offset(seq)
    off_idx = true_index(seq, i)
    chunk_idx = chunk(seq, i)
    bps = BioSequences.bits_per_symbol(Alphabet(seq)) % UInt
    first_carry_over = 5 - length(seq.data) % 5
    carry_over = 5 - (length(seq.data)-chunk_idx) % 5
    carried_over = 5-carry_over

    # println("\n$off_idx, $chunk_idx, $first_carry_over")
    # println("$((off_idx-1)*bps + first_carry_over), $(chunk_idx* 64)")
    # println("$carry_over, $carried_over")
    bound = (off_idx-1)*bps - chunk_idx*64 + first_carry_over
    # println(bound)
    cycled = (length(seq.data) - chunk_idx) % bps
    # println(cycled)

    if bound >= 0 && bound <= 64-bps
        # println("std: " * bin(seq.data[chunk_idx+1] << bound >> (64-bps)))
        return seq.data[chunk_idx+1] << bound >> (64-bps)
    else
        # println("except 1: " * bin((seq.data[chunk_idx]<<(64-carried_over)>>(64-carried_over-carry_over)) | (seq.data[chunk_idx+1] >> (64-carry_over))))
        # cycled_correction && println(cycled_correction)
        left_chunk_extract = seq.data[chunk_idx]<<(64-carried_over)>>(64-carried_over-carry_over)
        right_chunk_extract = seq.data[chunk_idx+1] >> (64-carry_over)
        return  left_chunk_extract | right_chunk_extract 
    end
end