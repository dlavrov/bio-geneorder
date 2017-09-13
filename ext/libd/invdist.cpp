#include "invdist.h"

int
calculate_offset ( Genome *g1, Genome *g2)
{
    intArray *genes2 = g2[0].pi;
    int num_genes = g2[0].len;
    int gA = g1[0].pi[0];

    for ( int i = 0; i < num_genes; i++ )
    {
        if ( genes2[i] == gA )
            return ( i );
        if ( genes2[i] == -gA )
            return ( i + num_genes );
    }

    return ( ERR_CONTENT );
}

int
num_breakpoints ( int *perm, int size )
{
    int i, b;
    int pA, pA1;
    int pB, pB1;

    b = 0;
    if ( perm[1] != 1 )
        b++;
    for ( i = 2; i < size - 1; i += 2 )
    {
        pA = perm[i];
        pA1 = ( pA == size ? 1 : pA + 1 );
        pB = perm[i + 1];
        pB1 = ( pB == size ? 1 : pB + 1 );
        if ( ( pB != pA1 ) && ( pA != pB1 ) )
            b++;
    }

    return ( b );
}

int
num_cycles ( int *perm, int size, distmem_t * distmem )
{
    int *done = distmem->done;
    int *greyEdges = distmem->greyEdges;
    int *invperm = distmem->done;
    int *cycle = distmem->labeled;
    
    int c = 0;
	int i, ind, j1, j2;
    int next;
    
#if 0
    /* This is now implicit */
    /* set black -edges */
    for ( i = 0; i < size; i += 2 )
    {
        blackEdges[i] = i + 1;
        blackEdges[i + 1] = i;
    }
#endif

    /* set grey -edges */
    for ( i = 0; i < size; i++ )
    {
        invperm[perm[i]] = i;
        greyEdges[i] = -1;
    }

    j1 = invperm[1];
    if ( j1 != 1 )
        greyEdges[0] = j1;

    for ( i = 1; i < size - 1; i += 2 )
    {
        ind = perm[i];
        if ( ind < perm[i + 1] )
        {
            j1 = invperm[ind - 1];
            j2 = invperm[ind + 2];
        }
        else
        {
            j1 = invperm[ind + 1];
            j2 = invperm[ind - 2];
        }
        if ( j1 != i - 1 )
        {
            greyEdges[i] = j1;
        }
        if ( j2 != i + 2 )
        {
            greyEdges[i + 1] = j2;
        }
    }

    j1 = invperm[size - 2];
    if ( j1 != size - 2 )
        greyEdges[size - 1] = j1;

    for ( i = 0; i < size; i++ )
        done[i] = 0;

    for ( i = 0; i < size; i++ )
    {
        if ( done[i] == 0 && greyEdges[i] != -1 )
        {
            cycle[i] = i;
            done[i] = 1;
            next = i;
            do
            {
                if ( next % 2 == 0 )
                    next++;
                else
                    next--;
                done[next] = 1;
                cycle[next] = i;
                next = greyEdges[next];
                done[next] = 1;
                cycle[next] = i;
            }
            while ( next != i );
            c++;
        }
    }

    return ( c );
}

void
connected_component ( int size, distmem_t * distmem, int *num_components )
{
    int i;
    int stack_ptr;
    int *greyEdges;
    int *cycle, *range, *cc;
    component_t *components;
    int *stack_root, *stack_range;
    int *parent, *next;
    int right, p;

    greyEdges = distmem->greyEdges;

    next = stack_root = distmem->stack;
    range = distmem->oriented;
    stack_range = cc = distmem->cc;
    parent = cycle = distmem->labeled;
    components = distmem->components;
    stack_ptr = -1;
    *num_components = 0;

    /* Use Linear algorithm to compute connected component */
    for ( i = 0; i < size; i++ )
    {
        if ( greyEdges[i] == -1 )
            continue;
        range[cycle[i]] = i;
    }

    for ( i = 0; i < size; i++ )
    {
        if ( greyEdges[i] == -1 )
            continue;           /*it is self loop,discard it */
        if ( parent[i] == i )
        {
            stack_ptr++;
            stack_root[stack_ptr] = i;
            stack_range[stack_ptr] = range[i];
        }
        else
        {                       /*check the top of stack for intersection */
            right = i;
            while ( stack_root[stack_ptr] > parent[i] )
            {
                /*union top to the i's connected component */
                parent[stack_root[stack_ptr]] = parent[i];
                if ( right < stack_range[stack_ptr] )
                    right = stack_range[stack_ptr]; /*extend the active range */
                stack_ptr--;
            }
            if ( stack_range[stack_ptr] < right )
                stack_range[stack_ptr] = right;
            if ( stack_range[stack_ptr] <= i )
            {
                /*the top connected-component is INACTIVE */
                components[*num_components].index = stack_root[stack_ptr];
                ( *num_components )++;
                stack_ptr--;
            }
        }
    }
    /*turn the forest to set of linked lists whose list head is index of
       component */
    for ( i = 0; i < size; i++ )
        next[i] = -1;
    for ( i = 0; i < size; i++ )
    {
        if ( greyEdges[i] == -1 )
            cc[i] = -1;
        else if ( i != parent[i] )
        {
            /* insert i between parent(i) and next of parent(i) */
            next[i] = next[parent[i]];
            next[parent[i]] = i;
        }
    }

    /*label each node with its root */
    for ( i = 0; i < *num_components; i++ )
    {
        p = components[i].index;
        while ( p != -1 )
        {
            cc[p] = i;
            p = next[p];
        }
    }

    return;
}

void
num_hurdles_and_fortress ( int *perm, int size,
                           int *num_hurdles, int *num_fortress,
                           distmem_t * distmem )
{
    int cIdx;
    int i, j;
    int *oriented, *cc, *labeled;
    component_t *components;
    int num_components;
    int num_oriented;
    int first_comp, last_comp, num_block;
    int num_superhurdles;
    int *greyEdges;

    /* By default, set number of hurdles and fortresses to 0 */
    *num_hurdles = 0;
    *num_fortress = 0;

    greyEdges = distmem->greyEdges;

    oriented = distmem->oriented;
    cc = distmem->cc;
    labeled = distmem->labeled;
    components = distmem->components;

    connected_component ( size, distmem, &num_components );


    if ( num_components == 0 )
    {
        return;
    }


    /* Calculate if each gray edge is oriented or unoriented */
    /* At the same time, label the connected component of a vertex as the
       index of its root */


    for ( i = 0; i < size; i++ )
    {
        j = greyEdges[i];
        if ( j == -1 )
        {
            oriented[i] = false;
        }
        else
        {
            if ( i < j )
            {
                if ( ( j - i ) % 2 != 0 )
                {
                    oriented[i] = false;
                    oriented[j] = false;
                }
                else
                {
                    oriented[i] = true;
                    oriented[j] = true;
                }
            }
        }
    }

    /* Look for any vertices that are "labeled".
       These are the roots of the cnnected components.
       Record them in the "components" array,
       and set num_components.
     */


    /* If a component contains an oriented vertex then it is oriented.
       Otherwise, the component is unoriented. */

    for ( i = 0; i < num_components; i++ )
    {
        components[i].oriented = false;
    }

    for ( i = 0; i < size; i++ )
    {
        if ( oriented[i] == 1 )
            components[cc[i]].oriented = true;
    }


    /* Count nonoriented components */

    num_oriented = 0;
    for ( i = 0; i < num_components; i++ )
    {
        if ( components[i].oriented == false )
        {
            num_oriented++;
        }
    }

    if ( num_oriented == 0 )
    {
        return;
    }

    for ( i = 0; i < num_components; i++ )
    {
        components[i].blocks = 0;
        components[i].hurdle = 0;
        components[i].left = -1;
        components[i].right = -1;
    }

    /* HURDLES 
       Hurdles are a subset of the nonoriented components. 
       There are two types of hurdles (in the KST-sense):
       "simple" and "superhurdle".
       First, we implicitly eliminate oriented components.
       Second, if a nonoriented component is one contiguous block of
       vertices, it is a hurdle.
       Third, if a hurdle "protects" the same non-hurdle on its left and
       right side, then it is a "superhurdle".
     */

    first_comp = -1;
    last_comp = -1;
    num_block = -1;
    for ( i = 0; i < size; i++ )
    {
        cIdx = cc[i];
        if ( cIdx != -1 )
        {
            if ( components[cIdx].oriented == false )
            {
                if ( cIdx != last_comp )
                {
                    if ( last_comp == -1 )
                    {
                        first_comp = cIdx;
                    }
                    else
                    {
                        components[last_comp].right = cIdx;
                        components[cIdx].left = last_comp;
                    }
                    last_comp = cIdx;
                    num_block++;
                    components[cIdx].blocks++;
                }
            }
        }
    }


    for ( i = 0; i < num_components; i++ )
    {
        if ( ( components[i].oriented == false )
             && ( components[i].blocks == 1 ) )
        {
            components[i].hurdle = HURDLE;
            ( *num_hurdles )++;
        }
    }
    if ( ( first_comp == last_comp )
         && ( components[first_comp].blocks == 2 ) )
    {
        components[first_comp].hurdle = ( HURDLE | GREATHURDLE );
        ( *num_hurdles )++;
    }

    if ( *num_hurdles < 3 )
        return;

    num_superhurdles = 0;
    for ( i = 0; i < num_components; i++ )
    {
        if ( components[i].hurdle )
        {
            if ( ( components[i].left == components[i].right ) &&
                 ( components[i].left != -1 ) )
            {
                if ( ( components[components[i].left].blocks == 2 ) &&
                     ( ( components[components[i].left].
                         hurdle & GREATHURDLE ) == 0 ) )
                {
                    components[i].hurdle |= SUPERHURDLE;
                    num_superhurdles++;
                }
                else
                {
                    return;
                }
            }
            else
            {
                return;
            }
        }
    }


    /* Set num_fortress if there are an odd number of hurdles,
       all of which are superhurdles. */
    if ( ( *num_hurdles == num_superhurdles )
         && ( num_superhurdles % 2 == 1 ) )
        *num_fortress = 1;

    return;
}

int
invdist_noncircular ( Genome * g1, Genome * g2, int offset )
{
    int i, twoi;
    int b, c;
    int g;
    int num_hurdles;
    int num_fortress;
    int reversal_dist;
    
    int num_genes = g1[0].len;
    int n = 2 * num_genes + 2;
    
    distmem_t * distmem = new distmem_t(n);

    int *perm1 = distmem->perm1;
    int *perm2 = distmem->perm2;
    int *perm = distmem->perm;

    for ( i = 0; i < num_genes; i++ )
    {
        g = 2 * g1[0].pi[i];
        twoi = 2 * i + 1;
        if ( g > 0 )
        {
            perm1[g - 1] = twoi;
            perm1[g] = twoi + 1;
        }
        else
        {
            perm1[-g] = twoi;
            perm1[-g - 1] = twoi + 1;
        }
    }
    
    for ( i = 0; i < num_genes; i++ )
    {
        if ( offset < num_genes )
            g = g2[0].pi[( offset + i ) % num_genes];
        else
        {
            g = -g2[0].pi[( offset - i ) % num_genes];
        }
        twoi = 2 * i + 1;
        if ( g > 0 )
        {
            perm2[twoi] = 2 * g - 1;
            perm2[twoi + 1] = 2 * g;
        }
        else
        {
            perm2[twoi] = -2 * g;
            perm2[twoi + 1] = -2 * g - 1;
        }
    }
    perm[0] = 0;
    for ( i = 1; i < n - 1; i++ )
    {
        perm[i] = perm1[perm2[i]];
    }
    perm[n - 1] = n - 1;

    b = num_breakpoints ( perm, n );
    c = num_cycles ( perm, n, distmem );
    num_hurdles_and_fortress ( perm, n, &num_hurdles, &num_fortress,
                               distmem );
    reversal_dist = b - c + num_hurdles + num_fortress;
    delete distmem;

    return ( reversal_dist );
}

int
invdist_circular ( Genome * g1, Genome * g2 )
{

    int offset;
	
    offset = calculate_offset ( g1, g2 );

    return ( invdist_noncircular ( g1, g2, offset ) );
}