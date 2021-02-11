### Commands used for openMS feature identification


## PeakPickerHiRes
## FeatureFinderCentroided
ls -1 ./*mzML |  parallel --plain -j 2 "./OpenMS-2.5.0/bin/PeakPickerHiRes -in {} -out ../1.peakpicker/{.}.ms1cent.mzML -algorithm:signal_to_noise 0 -algorithm:ms_levels 1 && ./OpenMS-2.5.0/bin/FeatureFinderCentroided -in ../1.peakpicker/{.}.ms1cent.mzML -out ../2.MS1features/{.}.featureXML"

##MapAlignerPoseClustering

./OpenMS-2.5.0/bin/MapAlignerPoseClustering -in `cat ./2.MS1features/filenames.txt` -out `cat ./3.AlignedFeat/filenames.txt` -algorithm:pairfinder:distance_RT:max_difference 300.0 -algorithm:pairfinder:distance_MZ:max_difference 20.0 -algorithm:pairfinder:distance_MZ:unit ppm

## FeatureLinkerUnlabeledQT
./OpenMS-2.5.0/bin/FeatureLinkerUnlabeledQT -in `cat ./3.AlignedFeat/filenames.txt` -out ./4.LinkedFeat/LinkedFeatures.consensusXML  -algorithm:distance_RT:max_difference 300.0 -algorithm:distance_MZ:max_difference 20.0 -algorithm:distance_MZ:unit ppm

## TextExporter
./OpenMS-2.5.0/bin/TextExporter -in ./4.LinkedFeat/LinkedFeatures.consensusXML -out ./4.LinkedFeat/ms1features.csv -quoting none  -id:peptides_only  -id:add_metavalues -1 -id:add_hit_metavalues -1 -id:add_protein_hit_metavalues -1 -consensus:sorting_method none

