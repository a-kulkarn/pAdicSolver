struct TropicalInterval
    min::Hecke.fmpq
    max::Union{Hecke.fmpq, AbstractFloat}
end

import Base: min, max, minimum, maximum, ==

# Access functions.
min(x::TropicalInterval) = x.min
infimum(::TropicalInterval) = min
minimum(::TropicalInterval) = min

max(x::TropicalInterval) = x.max
maximum(::TropicalInterval) = max
supremum(::TropicalInterval) = max


# Comparison.
(==)(x::TropicalInterval, y::TropicalInterval) = min(x) == min(y) && max(x) == max(y)

