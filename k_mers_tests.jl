using Pkg
Pkg.activate("k_mers/")

using Kmers
using Random

sample(l::N) where {N<:AbstractArray} = l[rand((1:length(l)))]

function tests()
    nucleotides = String["A", "T", "G", "C"]
    amino_acides = String["G", "A", "V", "L", "I", "T", "S", "M", "C", "P", "F", "Y", "W",
                          "H", "K", "R", "D", "E", "N", "Q"]

    testDNA31mer = join([sample(nucleotides) for _ in 1:31])
    testAA31mer = join([sample(amino_acides) for _ in 1:31])
    println("Test DNAKmer: $testDNA31mer")
    println("Test AAKmer: $testAA31mer")

    println("Kmers package: DNAKmer method")
    @time dnamer = DNAKmer{31}(testDNA31mer)
    println(Base.summarysize(dnamer))

    println("Kmers package: AAKmer method")
    @time aamer = AAKmer{31}(testAA31mer)
    println(Base.summarysize(aamer))
end

tests()