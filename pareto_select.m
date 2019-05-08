function [dual_obj,obj_tobecluster] = pareto_select(interval,lower,upper,scale_density,scale_distance,density_threshold,distance_threshold,dmax,d_list,unique_ID,directed_graph,density,similarity_density,similarity_distance)

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
        [H, iter, obj] = symnmf_newton(composite_similarity,2);
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