#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif

#include "distances.h"

Genome * structify(AV * pi) {
	int len;
	int i;
	SV ** elem;
	Genome * pi_perm;
	Genome perm;
	
	len = av_len(pi) +1;
	pi_perm = new Genome[len];
		
	for(i=0;i<len;i++){
		elem = av_fetch(pi,i,0);
		perm.pi = (intArray *)SvPV_nolen( *elem );
		perm.len = SvCUR( *elem ) / sizeof(intArray) -1;
		perm.circular = perm.pi[0];
		perm.pi++;
		
		pi_perm[i] = perm;
	}
	
	return pi_perm;
}

MODULE = Bio::GeneOrder::Distance		PACKAGE = Bio::GeneOrder::Distance

PROTOTYPES: ENABLE

intArray *
adjacencies_xs(pi,id)
	AV * pi
	AV * id
	CODE:
		Genome * pi_genome = structify(pi);
		Genome * id_genome = structify(id);
		
		std::vector<intArray> shared = _adjacencies(pi_genome,id_genome);
		
		delete pi_genome;
		delete id_genome;

		int size_RETVAL = shared.size();
		RETVAL = &shared[0];
	OUTPUT:
		RETVAL

int
breakpoints_xs(pi,id)
	AV * pi
	AV * id
	CODE:
		Genome * pi_genome = structify(pi);
		Genome * id_genome = structify(id);
		
		RETVAL = _breakpoints(pi_genome,id_genome);
		
		delete pi_genome;
		delete id_genome;
	OUTPUT:
		RETVAL
		
int
inversions_xs(pi,id)
	AV * pi
	AV * id
	CODE:
		Genome * pi_genome = structify(pi);
		Genome * id_genome = structify(id);
		
		RETVAL = _inversions(pi_genome,id_genome);
		
		delete pi_genome;
		delete id_genome;
	OUTPUT:
		RETVAL