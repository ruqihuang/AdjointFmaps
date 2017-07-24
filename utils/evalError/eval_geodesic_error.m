function error = eval_geodesic_error(S, match, gt_match)

	match = reshape(match, length(match), 1);
    gt_match = reshape(gt_match, size(match)); 
	pairs = [match gt_match];

	error = dijkstra_pairs(S, pairs);
	error = error/S.sqrt_area;

end