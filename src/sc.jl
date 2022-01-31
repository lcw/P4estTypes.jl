module SC

using ..MPI
using ..P4est
initialized() = sc_package_id()[] >= 0
finalize() = sc_finalize_noabort()
module LP

using ..P4est
abstract type LogPriority end
struct Always <: LogPriority end
struct Debug <: LogPriority end
struct Default <: LogPriority end
struct Error <: LogPriority end
struct Essential <: LogPriority end
struct Info <: LogPriority end
struct Production <: LogPriority end
struct Silent <: LogPriority end
struct Statistics <: LogPriority end
struct Trace <: LogPriority end
struct Verbose <: LogPriority end

Base.convert(::Type{Cint}, ::Always) = Cint(SC_LP_ALWAYS)
Base.convert(::Type{Cint}, ::Debug) = Cint(SC_LP_DEBUG)
Base.convert(::Type{Cint}, ::Default) = Cint(SC_LP_DEFAULT)
Base.convert(::Type{Cint}, ::Error) = Cint(SC_LP_ERROR)
Base.convert(::Type{Cint}, ::Essential) = Cint(SC_LP_ESSENTIAL)
Base.convert(::Type{Cint}, ::Info) = Cint(SC_LP_INFO)
Base.convert(::Type{Cint}, ::Production) = Cint(SC_LP_PRODUCTION)
Base.convert(::Type{Cint}, ::Silent) = Cint(SC_LP_SILENT)
Base.convert(::Type{Cint}, ::Statistics) = Cint(SC_LP_STATISTICS)
Base.convert(::Type{Cint}, ::Trace) = Cint(SC_LP_TRACE)
Base.convert(::Type{Cint}, ::Verbose) = Cint(SC_LP_VERBOSE)

struct SCArray{T} <: AbstractArray{T,1}
    pointer::Ptr{sc_array}
end
Base.IndexStyle(::SCArray) = IndexLinear()
function Base.getindex(a::SCArray{T}, i::Int) where {T}
    return GC.@preserve a unsafe_load(T(a.pointer.array), i)
end
Base.size(a::SCArray) = (a.pointer.elem_count,)

end # module

setverbosity(lp::SC.LP.LogPriority) = sc_package_set_verbosity(sc_package_id()[], lp)

end # module
