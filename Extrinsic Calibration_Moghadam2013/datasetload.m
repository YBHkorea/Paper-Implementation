function [images, MLimages, pc_origin, objPCSet, cameraParams, GT_RT] = datasetload(dataset_num)
if dataset_num == 6
    dataset = "dataset" + dataset_num + ".mat";
elseif dataset_num == 7
    Path = Path + "/C304_4_nopattern";
end

load(dataset);

end