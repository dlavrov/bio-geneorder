#include "distances.h"
#include "invdist.h"

std::vector<intArray> _adjacencies(Genome * pi, Genome * id){

	int num_pi;
	int num_id;
	int i,j,k,l;
	
	num_pi = sizeof(pi)/sizeof(Genome *);
	num_id = sizeof(id)/sizeof(Genome *);
	
	std::vector<intArray> shared_bounds;
	
	// For each permutation in the identity genome
	for(i=0;i<num_id;i++){
		
		// Go through each permutation in the comparison genome
		for(j=0;j<num_pi;j++){
			int find;
			
			// And go through each gene boundary in the comparison permutation
			for(k=0;k<pi[j].len-1;k++){
				find = 0;
				
				// These are the genes in the boundary
				int pi1 = pi[j].pi[k];
				int pi2 = pi[j].pi[k+1];
				
				// Check if each boundary is in the identiy permutation
				for(l=0;l<id[i].len-1;l++){
					int id1 = id[i].pi[l];
					int id2 = id[i].pi[l+1];
					
					if( (pi1 == id1 && pi2 == id2) || (-pi1 == id2 && -pi2 == id1) ){
						find++;
					}
				}
				
				if(find){
					shared_bounds.push_back(pi1);
					shared_bounds.push_back(pi2);
				}
			}

			// If the comparison permutation is circular then check the ends
			if(pi[j].circular){
				find = 0;
				
				// These are the genes in the boundary
				int pi1 = pi[j].pi[pi[j].len-1];
				int pi2 = pi[j].pi[0];
				
				// Check if the end boundary is in the identity permutation
				for(l=0;l<id[i].len-1;l++){
					int id1 = id[i].pi[l];
					int id2 = id[i].pi[l+1];
					
					if( (pi1 == id1 && pi2 == id2) || (-pi1 == id2 && -pi2 == id1) ){
						find++;
					}
				}
				
				// Check the ends of the identity if it is ciruclar
				if(id[i].circular){
					int id1 = id[i].pi[id[i].len-1];
					int id2 = id[i].pi[0];
					
					if( (pi1 == id1 && pi2 == id2) || (-pi1 == id2 && -pi2 == id1) ){
						find++;
					}
				}
				
				if(find){
					shared_bounds.push_back(pi1);
					shared_bounds.push_back(pi2);
				}
			}
			
		}
	}

	return shared_bounds;
}

int _breakpoints(Genome * pi, Genome * id){

	int num_pi;
	int num_id;
	int i,j,k,l,bA,bB,bID,b;
	bA = bB = b = bID = 0;
	
	num_pi = sizeof(pi)/sizeof(Genome *);
	num_id = sizeof(id)/sizeof(Genome *);
	
	// For each permutation in the identity genome
	for(i=0;i<num_id;i++){
		bID += bB;
		
		// Go through each permutation in the comparison genome
		for(j=0;j<num_pi;j++){
			int find;
			
			// And go through each gene boundary in the comparison permutation
			for(k=0;k<pi[j].len-1;k++){
				find = 0;
				bA++;
				
				// These are the genes in the boundary
				int pi1 = pi[j].pi[k];
				int pi2 = pi[j].pi[k+1];
				
				// Check if each boundary is in the identiy permutation
				bB = 0;
				for(l=0;l<id[i].len-1;l++){
					int id1 = id[i].pi[l];
					int id2 = id[i].pi[l+1];
					bB++;
					
					if( (pi1 == id1 && pi2 == id2) || (-pi1 == id2 && -pi2 == id1) ){
						find++;
					}
				}
				
				if(find){
					b++;
				}
			}

			// If the comparison permutation is circular then check the ends
			if(pi[j].circular){
				find = 0;
				bA++;
				bB++;
				
				// These are the genes in the boundary
				int pi1 = pi[j].pi[pi[j].len-1];
				int pi2 = pi[j].pi[0];
				
				// Check if the end boundary is in the identity permutation
				for(l=0;l<id[i].len-1;l++){
					int id1 = id[i].pi[l];
					int id2 = id[i].pi[l+1];
					
					if( (pi1 == id1 && pi2 == id2) || (-pi1 == id2 && -pi2 == id1) ){
						find++;
					}
				}
				
				// Check the ends of the identity if it is ciruclar
				if(id[i].circular){
					int id1 = id[i].pi[id[i].len-1];
					int id2 = id[i].pi[0];
					
					if( (pi1 == id1 && pi2 == id2) || (-pi1 == id2 && -pi2 == id1) ){
						find++;
					}
				}
				
				if(find){
					b++;
				}
			}
			
		}
	}
	
	// The number of breakpoints is the difference between the number of adjacencies
	// and the size of the largest genome;
	b = bA > bB ? bA - b : bB - b;

	return b;
}

int _inversions(Genome * pi, Genome * id){

	int i;
	int inversions;
	
	int num_pi = sizeof(pi)/sizeof(Genome *);
	int num_id = sizeof(id)/sizeof(Genome *);
	
	if(duplicates(pi) || duplicates(id)){
		return ERR_DUPLICATES;
	}else if(unequal_content(pi,id)){
		return ERR_CONTENT;
	}else{
	
		if(num_pi == num_id == 1){
		
			if(pi[0].circular){
			
				inversions = invdist_circular(id,pi) + 1 - id[0].circular;
			
			}else if(id[0].circular){
			
				inversions = invdist_circular(pi,id) + 1;
				
			}else{
			
				inversions = invdist_noncircular(pi,id,0);
			}
			
		}else{
		
			for(i=0;i<num_pi || i<num_id;i++){
				if(pi[i].circular || id[i].circular){
					return _DCJ(pi,id);
				}
			}
			
			//inversions = mcdist_noncircular(pi,id);
			
		}	
	}
	
	return inversions;
}


int _DCJ(Genome * pi, Genome * id){

	int i;
	int dcj;
	
	int num_pi = sizeof(pi)/sizeof(Genome *);
	int num_id = sizeof(id)/sizeof(Genome *);
	
	if(duplicates(pi) || duplicates(id)){
		return ERR_DUPLICATES;
	}else if(unequal_content(pi,id)){
		return ERR_CONTENT;
	}else{
	
		return ERR_NOTIMPL;
	}
	
	
	return dcj;
}

bool duplicates(Genome * pi){
	int num_pi;
	int i,j;
	
	num_pi = sizeof(pi)/sizeof(Genome *);
	
	std::vector<intArray> genes;
	
	// For each permutation in the genome
	for(i=0;i<num_pi;i++){
			
			// Go through each gene
			for(j=0;j<pi[i].len;j++){
				
				genes.push_back(abs(pi[i].pi[j]));
				
			}
	}
	
	std::sort( genes.begin(), genes.end() );
	
	for(i=0;i<genes.size()-1;i++){
		if(genes[i] == genes[i+1]){
			return 1;
		}
	}
	return 0;
}

bool unequal_content(Genome * pi, Genome * id){
	int num_pi;
	int num_id;
	int i,j;
	
	num_pi = sizeof(pi)/sizeof(Genome *);
	num_id = sizeof(id)/sizeof(Genome *);
	
	std::vector<intArray> genesA;
	std::vector<intArray> genesB;
	
	// For each permutation in the genome
	for(i=0;i<num_pi;i++){
			
			// Go through each gene
			for(j=0;j<pi[i].len;j++){
				
				genesA.push_back(abs(pi[i].pi[j]));
				
			}
	}
	
	// For each permutation in the genome
	for(i=0;i<num_id;i++){
			
			// Go through each gene
			for(j=0;j<id[i].len;j++){
				
				genesB.push_back(abs(id[i].pi[j]));
				
			}
	}
	
	std::sort( genesA.begin(), genesA.end() );
	std::sort( genesB.begin(), genesB.end() );
	
	if( genesA.size() != genesB.size() ){
		return 1;
	}
	
	for(i=0;i<genesA.size();i++){
		if(genesA[i] != genesB[i]){
			return 1;
		}
	}
	
	return 0;
}