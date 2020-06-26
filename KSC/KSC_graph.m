function W =  KSC_graph(feature,para)

    Z = unifiedcluster(feature,para.alphak,para.betak);
    W = BuildAdjacency(thrC(Z,para.rho));
    
end
