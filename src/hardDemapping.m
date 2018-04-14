
function [hardDemapperBits,hardDemapperSymbol] = hardDemapping(sym_labeling,symbol_set,rxSig)
	% input: the symbols and its labeling
	%        received signal

	% smallest euclidian distance
	%[m, l] = get_mapping_and_labeling(nb,labeling,mapPamQam);

	l_max=length(symbol_set);
	nb =log2(l_max);
	hardDecision = zeros(1,length(rxSig));

	distance = zeros(l_max,length(rxSig));
	for n = 1:l_max
		distance(n,:) = abs(rxSig-symbol_set(n));
	end
	[~,index_min_distance] = min(distance);

	for n = 1:l_max
		hardDecision(index_min_distance==n) = sym_labeling(n);
	end

	% same result but less sufficient:
	%for a = 1:length(rxSig)
	%	dec=l(1);
	%	dist = realmax;
	%	for n = 1:l_max
	%	   distneu = abs(rxSig(a)-l(n));
	%	   if ( distneu < dist)
	%	       dec=m(n);
	%	       dist=distneu;
	%	   end
	%	end
	%	hardDecision(a)=dec;
	%end
	hardDemapperSymbol = symbol_set(index_min_distance);
	hardDemapperBits = de2bi(hardDecision,nb);

end
