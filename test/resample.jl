using GeoData, Test, ArchGDAL, GeoFormatTypes, NCDatasets, Dates
using GeoData: resample

include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

@testset "resample" begin
    raster_path = maybedownload("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

    output_res = 0.0027
    output_crs = EPSG(4326)
    resample_method = "near"

    ## Resample cea.tif manually with ArchGDAL
    wkt = convert(String, convert(WellKnownText, output_crs))
    AG_output = ArchGDAL.read(raster_path) do dataset
        ArchGDAL.gdalwarp([dataset], ["-t_srs", "$(wkt)",
                                "-tr", "$(output_res)", "$(output_res)",
                                "-r", "$(resample_method)"]) do warped
            ArchGDAL.read(ArchGDAL.getband(warped, 1))
        end
    end

    ## Resample cea.tif using resample
    cea = GeoArray(GDALarray(raster_path))
    GD_output = resample(cea, output_res, crs = output_crs, method = resample_method)

    ## Compare the two
    @test AG_output == GD_output.data[:, :, 1]
    @test abs(step(dims(GD_output, Lat))) ≈ abs(step(dims(GD_output, Lon))) ≈ output_res

    @testset "snapped size and dim index match" begin
        snaptarget = GD_output
        snapped = resample(cea, snaptarget)
        @test size(snapped) == size(snaptarget)
        @test isapprox(index(snaptarget, Lat), index(snapped, Lat))
        @test isapprox(index(snaptarget, Lon), index(snapped, Lon))
    end

    @testset "snap resample with extra dimension" begin
        ceati = cat(cea, cea, cea, cea; dims=Ti(Date(2001, 1):Month(1):Date(2001, 4)))
        snaptarget = GD_output
        snapped = resample(ceati, snaptarget)
        @test size(snapped) == (123, 103, 1, 4)
        @test isapprox(index(snaptarget, Lat), index(snapped, Lat))
        @test isapprox(index(snaptarget, Lon), index(snapped, Lon))
    end

    @testset "snap stack" begin
        stack = GeoStack((a=cea, b=cea))
        snappedstack = resample(stack, snaptarget)
        @test snappedstack isa GeoStack
        @test size(snappedstack[:a]) == size(snappedstack[:b]) == (123, 103, 1)
        @test isapprox(index(snaptarget, Lat), index(snappedstack[:a], Lat))
        @test isapprox(index(snaptarget, Lon), index(snappedstack[:b], Lon))
    end

end
