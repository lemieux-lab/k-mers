using Kmers, Kmers.BioSequences

struct AAAlphabet <: Alphabet end

BioSequences.BitsPerSymbol(::AAAlphabet) = BioSequences.BitsPerSymbol{5}()

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
chunk(k::Km, i::Int) where {Km<:Kmer} = true_index(k, i)÷(12+carry(true_index(k, i)÷12))

@inline function BioSequences.extract_encoded_element(seq::Kmer{<:AAAlphabet}, i::Integer)
    off = offset(seq)
    off_idx = true_index(seq, i)
    chunk_idx = chunk(seq, i)
    bps = BioSequences.bits_per_symbol(Alphabet(seq)) % UInt
    prev_carry_over = length(seq.data) % 5
    carry_over = (5-length(seq.data)+chunk_idx) % 5
    carried_over = 5-carry_over

    # println("\n$off_idx, $chunk_idx, $prev_carry_over")
    # println("$(off_idx*bps), $(chunk_idx* 64 + prev_carry_over)")
    # println((off_idx*bps - chunk_idx*64 + prev_carry_over) < 64)
    # println(bin(seq.data[chunk_idx+1] << (off_idx*bps - chunk_idx * 64 -prev_carry_over) >> (64-bps)))
    if (off_idx*bps >= chunk_idx*64 + prev_carry_over) #&& ((off_idx*bps - chunk_idx*64 + prev_carry_over) < 64)
        return seq.data[chunk_idx+1] << (off_idx*bps - chunk_idx * 64 - prev_carry_over) >> (64-bps)
    else
        # println(bin(seq.data[chunk_idx]<<(64-carried_over)>>(64-carried_over-carry_over)))
        return (seq.data[chunk_idx]<<(64-carried_over)>>(64-carried_over-carry_over)) | (seq.data[chunk_idx+1] >> (64-carry_over))
    end
end