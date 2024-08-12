function resPath = HyAs1(sPos, gPos, obsInfo, N)

resolution = 0.05;
K = 1/resolution;
map_width = 60;
map_length = 8;
costVal = zeros(map_length*K,map_width*K);
costVal(1,:) = 1;costVal(end,:) = 1;
costVal(:,1) = 1;costVal(:,end) = 1;
costVal(1:1*K,30*K:60*K) = 1;
costVal(4*K:8*K,30*K:60*K) = 1;
% costVal(2.2*K:4*K,40*K:50*K) = 1;
obs_x = obsInfo(1) + 30;
obs_y = 4 - obsInfo(2);
obs_w = obsInfo(3)+1.2;
for i = 1:obs_w*K
    for j = 1:obs_w*K
        costVal(round(obs_y*K+i-obs_w*K/2),round(obs_x*K+j-obs_w*K/2)) = 1;
    end
end
map = binaryOccupancyMap(costVal,K);
% show(map);
ss = stateSpaceSE2;
ss.StateBounds = [map.XWorldLimits;map.YWorldLimits;[-pi pi]];
sv = validatorOccupancyMap(ss);
sv.Map = map;
planner = plannerHybridAStar(sv, ...
                             MinTurningRadius=2, ...
                             MotionPrimitiveLength=0.6);
sPos = sPos' + [30 4 0];
gPos = gPos' + [30 4 0];
refpath = plan(planner,sPos,gPos,SearchMode='exhaustive');
figure(1)
show(planner)

tmp_s = size(refpath.States,1);
if(tmp_s >= N)
    resPath = refpath.States(1:N,:);
else
    resPath = refpath.States;
    resPath = [resPath;repmat(resPath(end,:),[N-tmp_s 1])];
end
resPath = resPath - [30 4 0];
% show(planner)
% pause(0.5);
% close(gcf)


end