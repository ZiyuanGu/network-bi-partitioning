function [dual_obj,obj_tobecluster] = pareto_select(interval,lower,upper,scale_density,scale_distance,density_threshold,distance_threshold,dmax,d_list,unique_ID,directed_graph,density,similarity_density,similarity_distance)

% The function has several input arguments:
% interval -- the number of intervals between 0 and 1, which is also the increment of the weight coefficient
% lower and upper -- the lower and upper bounds of the weight coefficient, which are 0 and 1
% scale_density and scale_distance -- the scaling parameters for calculating the density and the distance similarity measures
% density_threshold and distance_threshold -- the thresholds for calculating the density and the distance similarity measures
% dmax -- the maximum shortest path from the source link to any other link in the network
% d_list -- the shortest paths from the source link to all the other links in the network
% unique_ID -- the unique external identifier for each link in the network
% directed_graph -- a matrix where (in our case) the 1st column is the unique external identifier (unique_ID), the 2nd, 3rd and 4th columns are the internal identifiers for each link and its two end nodes in AIMSUN
% density -- the density data
% similarity_density and similarity_distance -- two column vectors of the density and the distance similarity measures for each link in the network
% 
% The function has two output matrices:
% dual_obj -- a matrix where the 1st column is the incremental weight coefficient, the 2nd and 3rd (4th and 5th) columns are the average density (distance) simialrity measures for the two clusters
% obj_tobecluster -- a matrix where the 1st and 2nd columns are the average density and distance similarity measures for the RZ, the 3rd and 4th columns are their respective changes compared with the previous values, 
% the 5th column is the overall change which essentially gives the significant solutions from the Pareto front

dual_obj = [];

for w = lower:((upper-lower)/interval):upper
    similarity_index = zeros(945,1);
    assign_list = zeros(945,1);
    composite_similarity = zeros(945);
    w_density = w;
    w_distance = 1-w_density;
    
    for i = 1:945
        index = find(unique_ID(:,1) == directed_graph(i,1));
        
        if density(i,17) > density_threshold && d_list(index,1) < distance_threshold
            similarity_index(i,1) = 1;
            
        elseif density(i,17) <= density_threshold && d_list(index,1) < distance_threshold
            similarity_index(i,1) = w_distance+w_density*(density(i,17)/density_threshold).^scale_density;
            
        elseif density(i,17) > density_threshold && d_list(index,1) >= distance_threshold
            similarity_index(i,1) = w_density+w_distance*((dmax-d_list(index,1))/(dmax-distance_threshold)).^scale_distance;
            
        else
            similarity_index(i,1) = w_density*(density(i,17)/density_threshold).^scale_density+w_distance*((dmax-d_list(index,1))/(dmax-distance_threshold)).^scale_distance;
        end
    end
    
    for i = 1:945
        for j = 1:945
            composite_similarity(i,j) = 1-abs(similarity_index(i,1)-similarity_index(j,1));
        end
    end
    
    iterations_H = cell(10,1);
    iterations_iter = zeros(10,1);
    iterations_obj = zeros(10,1);
    
    for i =1:10
        [H,iter,obj] = symnmf_newton(composite_similarity,2);
        iterations_H{i,1} = H;
        iterations_iter(i,1) = iter;
        iterations_obj(i,1) = obj;
    end
    
    H = iterations_H{find(iterations_obj(:,1) == min(iterations_obj(:,1))),1};
    
    for i = 1:945
        [~,I] = max(H(i,:),[],2);
        assign_list(i,1) = I;
    end
    
    cluster_density_1 = [];
    cluster_density_2 = [];
    cluster_distance_1 = [];
    cluster_distance_2 = [];
    
    for i = 1:945
        if assign_list(i,1) == 1
            cluster_density_1 = [cluster_density_1;similarity_density(i,1)];
            cluster_distance_1 = [cluster_distance_1;similarity_distance(i,1)];
            
        else
            cluster_density_2 = [cluster_density_2;similarity_density(i,1)];
            cluster_distance_2 = [cluster_distance_2;similarity_distance(i,1)];
        end
    end
    
    dual_obj = [dual_obj;w,mean(cluster_density_1),mean(cluster_density_2),mean(cluster_distance_1),mean(cluster_distance_2)];
end

obj_tobecluster = zeros(size(dual_obj,1),5);

for i = 1:size(obj_tobecluster,1)
    obj_tobecluster(i,1:2) = [max(dual_obj(i,2),dual_obj(i,3)),max(dual_obj(i,4),dual_obj(i,5))];
end

for i = 2:size(obj_tobecluster,1)
    obj_tobecluster(i,3) = 100*(obj_tobecluster(i,1)-obj_tobecluster(i-1,1))/obj_tobecluster(i-1,1);
    obj_tobecluster(i,4) = 100*(obj_tobecluster(i,2)-obj_tobecluster(i-1,2))/obj_tobecluster(i-1,2);
    obj_tobecluster(i,5) = obj_tobecluster(i,3)+obj_tobecluster(i,4);
end
end