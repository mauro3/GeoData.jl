export resample

# Dispatch trait to process non-gdal dimensions
struct ExtraDims end
struct GDALdims end

dimsavailable(A) =
    length(otherdims(A, (Lat(), Lon(), Band()))) > 0 ? ExtraDims() : GDALdims()

"""
	resample(x, resolution::Number; crs::GeoFormat=crs(A), method::String="near")
	resample(x, snap::AbstractGeoArray; method::String="near")

`resample` uses `ArchGDAL.gdalwarp` to resample any `AbstractGeoArray`.
`Dimension`s other then `Lat`, `Lon` and `Band` are sliced and processed 
separately, then rejoined into a single array. When `x` is and `AbstrackGeoStack` 
or `AbstractGeoSeries` all child data is resampled.

## Arguments
- `x`: An `AbstractGeoArray`, `AbstractGeoStack` or `AbstractGeoSeries` to resample.
- `resolution`: A `Number` specifying the resolution for the output. 
  If the keyword argument `crs` (described below) is specified, `resolution` must be in units of the `crs`.
- `snap`: an `AbstractGeoArray` whos resolution, crs and bounds will be snapped to. 
  For best results it should roughly cover the same extent, or a subset of `A`.

## Keyword Arguments
- `crs`: A `GeoFormatTypes.GeoFormat` specifying an output crs 
  (`A` with be reprojected to `crs` in addition to being resampled). Defaults to `crs(A)`
- `method`: A `String` specifying the method to use for resampling. Defaults to `"near"` 
  (nearest neighbor resampling). See [resampling method](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r) 
  in the gdalwarp docs for a complete list of possible values.

"""
resample

resample(s::Union{AbstractGeoSeries,AbstractGeoStack}, args...; kw...) =
    map(x -> resample(x, args...; kw...), s)

resample(A::AbstractGeoArray, args...; kw...) =
    _resample(dimsavailable(A), A, args...; kw...)

# Slice `A` using the `dimwise` generator form DD, and 
# resample each slice. Then `cat` them back together.
function _resample(::ExtraDims, A::AbstractGeoArray, args...; kw...)
    othdims = otherdims(A, (Lat(), Lon(), Band()))
    slices = map(DD.dimwise_generators(othdims)) do d
        _resample(GDALdims(), A[d...], args...; kw...)
    end
    cat(slices...; dims=othdims)
end

function _resample(::GDALdims, A::AbstractGeoArray, resolution::Number;
				  crs::GeoFormat=crs(A),
                  method::String="near")
    wkt = convert(String, convert(WellKnownText, crs))
    flags = ["-t_srs", "$(wkt)",
             "-tr", "$(resolution)", "$(resolution)",
             "-r", "$(method)"]
    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], flags) do warped
            GeoArray(warped)
        end
    end
end

function _resample(::GDALdims, A::AbstractGeoArray, snap::AbstractGeoArray; method::String="near")
    wkt = convert(String, convert(WellKnownText, crs(snap)))
    latres, lonres = map(abs ∘ step, span(snap, (Lat(), Lon())))
    (latmin, latmax), (lonmin, lonmax) = bounds(snap, (Lat(), Lon()))
    flags = ["-t_srs", "$(wkt)",
             "-tr", "$(latres)", "$(lonres)",
             "-te", "$(lonmin)", "$(latmin)", "$(lonmax)", "$(latmax)",
             "-r", "$(method)"]
    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], flags) do warped
            GeoArray(warped)
        end
    end
end
