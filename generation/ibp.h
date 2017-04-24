**
* @file ibp.h
*
* Routines for generating IBP relations.
*
#ifndef `IBPHFILE'
#define IBPHFILE
*--#[ decl :

**
* Internals, Externals, Propagators, Replacements: the input of PrepareIBP().
*
#define Internals "UNDEF"
#define Externals "UNDEF"
#define Propagators ""
#define Replacements ""

**
* The output of PrepareIBP().
*
#define NumberOfLoops "0"
#define NumberOfPropagators "0"

**
* The output of GenerateIBP().
*
#define NumberOfRelations "0"

**
* Internally used preprocessor variables.
*
#define NormalizedInternals ""
#define NormalizedExternals ""
#define NormalizedPropagators ""
#define NormalizedReplacements ""
#define IBPYOperator ""
#define IdentifySPs ""
#define IdentifyYOps ""

*--#] decl :
*--#[ proc :
*--#[ PrepareIBP :

**
* Prepares for the following IBP relation generations. It uses preprocessor
* variables Internals, Externals and Propagators as the input. The argument
* "Y" is used for declaring symbols Yi meaning Yi^n decreases the i-th index
* by n.
*
* Example:
*   V k1,k2,p;
*   #define Internals   "k1"
*   #define Externals   "p"
*   #define Propagators "k1.k1,k2.k2"
*   #define Replacements "
*     id k2 = k1 + p;
*   "
*   #call PrepareIBP(Y)
*
* This procedure contains
*   - Declaration statements.
*
#procedure PrepareIBP(Y)
  #define prefix "IBPPrivate"
* Internally used functions.
  #define R      "`prefix'R"
  #define E      "`prefix'E"
  #define G      "`prefix'G"
* Internally used $-variables.
  #define tmp    "`prefix'tmp1"
  #define tmp1   "`prefix'tmp1"
  #define tmp2   "`prefix'tmp2"
  #define FF     "`prefix'FF"
  #define P      "`prefix'P"
  #define found  "`prefix'found"
  #define ss     "`prefix'ss"
  #define coeff  "`prefix'coeff"
  #define rhs    "`prefix'rhs"
  #define res    "`prefix'res"

  #if "`Internals'" == "UNDEF"
    #message Error: Internals must be defined
    #terminate
  #endif

  #if "`Externals'" == "UNDEF"
    #message Error: Externals must be defined
    #terminate
  #endif

  V `Internals',`Externals';  * Must have valid names.
  CF `R',`E',`G';

* First, normalize the input. We allow the input to contain line breaks, because
* it is convenient especially for Propagators, which can easily become long.
* Line breaks cause troubles if they are used as #do s={`input'} (#do discards
* the text after the first line break), so we need to remove such line breaks
* from the input.

* Normalize Internals.
  #$`tmp1' = `R'(`Internals');
  #inside $`tmp1'
    if (match(`R'(?a$`tmp2')));
    endif;
  #endinside
  #define Internals "`$`tmp2''"

* Normalize Externals.
  #$`tmp1' = `R'(`Externals');
  #inside $`tmp1'
    if (match(`R'(?a$`tmp2')));
    endif;
  #endinside
  #define Externals "`$`tmp2''"

* Normalize Propagators.
  #$`tmp1' = `R'(`Propagators');
  #inside $`tmp1'
    if (match(`R'(?a$`tmp2')));
    endif;
  #endinside
  #define Propagators "`$`tmp2''"

* Check if Replacements is valid.
  #$`tmp' = 1;
  #inside $`tmp1'
    `Replacements';
  #endinside

* Store the input for later use.
  #redefine NormalizedInternals "`Internals'"
  #redefine NormalizedExternals "`Externals'"
  #redefine NormalizedPropagators "`Propagators'"
  #redefine NormalizedReplacements "`Replacements'"
  #redefine IBPYOperator "`Y'"

* Print the input.
  #ifdef `V'
    #message Internals = {`Internals'}
    #message Externals = {`Externals'}
    #message Propagators = {`Propagators'}
    #write "~~~Replacements = {`Replacements'}"
  #endif

* Check Internals.
  #define NI "0"
  #do k={`Internals',}
    #ifdef `k'
      #redefine NI "{`NI'+1}"
      #define Internal`NI' "`k'"
    #endif
  #enddo
  #if `NI' == 0
    #message Error: Internals must be defined
    #terminate
  #endif

* Check Externals. Externals can be empty.
  #define NE "0"
  #do k={`Externals',}
    #ifdef `k'
      #redefine NE "{`NE'+1}"
      #define External`NE' "`k'"
    #endif
  #enddo

* Check Propagators.
  #define NP "0"
  #do P={`Propagators',}
    #ifdef `P'
      #redefine NP "{`NP'+1}"
      #define Propagator`NP' "`P'"
    #endif
  #enddo
  #if `NP' == 0
    #message Error: Propagators must be defined
    #terminate
  #endif

* Check if all scalar products can be expressed in terms of the propagators.
  #define NS "{`NI'*(`NI'+1)/2 + `NI'*`NE'}"
  #if `NS' < `NP'
    #message Error: linearly dependent propagators
    #terminate
  #else if `NS' > `NP'
    #message Error: not enough propagators
    #terminate
  #endif

* Construct a list of scalar products containing loop momenta.
  #redefine NS "0"
  #define SPs ""
  #do i=1,`NI'
    #do j=`i',`NI'
      #redefine NS "{`NS'+1}"
      #$`tmp' = `Internal`i''.`Internal`j'';
      #redefine SPs "'SPs'`$`tmp'',"
      #define SPs`NS' "`$`tmp''"
    #enddo
  #enddo
  #do i=1,`NI'
    #do j=1,`NE'
      #redefine NS "{`NS'+1}"
      #$`tmp' = `Internal`i''.`External`j'';
      #redefine SPs "'SPs'`$`tmp'',"
      #define SP`NS' "`$`tmp''"
    #enddo
  #enddo

* Solve the linear system of equations Yi = P(i) for the scalar products,
* where Yi is an operator decreasing the i-th index by one, and P(i) is
* the i-th propagator. Each propagator is a linear combination of the scalar
* products and other constants, therefore the scalar products must eventually
* be expressed as linear combinations of operators Yi and constants.

  S `Y'1,...,`Y'`NP';

* $FF stores the solutions such that the result of "s" is the coefficient
* of E(s). G(s) is used for holding unsolved variables.
  #$`FF' =
    #do s={`SPs',}
      #ifdef `s'
        + `E'(`s') * `G'(`s')
      #endif
    #enddo
  ;
  #do i=1,`NP'
*   Consider one of the propagators.
    #$`P' = `Propagator`i'';
    #inside $`P'
      repeat;
        `NormalizedReplacements';
      endrepeat;
    #endinside
    #do s={`SPs',}
      #ifdef `s'
*       Effectively, id s = $FF[E(s)] in $P.
        #$`tmp' = $`FF';
        #inside $`tmp'
          id `E'(`s') = 1;
          id `E'(?a) = 0;
        #endinside
        #inside $`P'
          id `s' = `$`tmp'';
        #endinside
      #endif
    #enddo
    #inside $`P'
      repeat;
        `NormalizedReplacements';
      endrepeat;
    #endinside
*   Pick up one of unsolved scalar products.
    #$`found' = 0;
    #inside $`P'
      #do s={`SPs',}
        #ifdef `s'
          if ($`found' == 0);
            if (match(`G'(`s')));
              $`found' = 1;
              $`ss' = `s';
            endif;
          endif;
        #endif
      #enddo
    #endinside
    #if `$`found'' == 0
      #message Error: lsolve failed (`Y'`i' == `$`P'')
      #terminate
    #endif
*   Solve Yi = $P with respect to G($ss).
    #$`coeff' = $`P';
    #inside $`coeff'
      id ifmatch -> keep`i', `G'($`ss') = 1;
      discard;
label keep`i';
    #endinside
    #$`rhs' = $`P';
    #inside $`rhs'
      id `G'($`ss') = 0;
    #endinside
*   XXX: $coeff is assumed to be simple and no denom_ in the result.
    #$`res' = (`Y'`i' - $`rhs') / $`coeff';
*   Store the result.
    #inside $`FF'
      id `G'($`ss') = $`res';
    #endinside
  #enddo
* Build id statements for replacing the scalar products.
  #redefine IdentifySPs ""
  #do s={`SPs',}
    #ifdef `s'
*     Add a rule, effectively, id s = `$FF[E(s)]'.
      #$`tmp' = $`FF';
      #inside $`tmp'
        id `E'(`s') = 1;
        id `E'(?a) = 0;
      #endinside
      #redefine IdentifySPs "`IdentifySPs'
  id `s'=`$`tmp'';"
    #endif
  #enddo
  #ifdef `V'
    #write "~~~IdentifySPs = {`IdentifySPs'
}"
  #endif
* Redefine variables.
  #redefine NumberOfLoops "`NI'"
  #redefine NumberOfPropagators "`NP'"
#endprocedure

*--#] PrepareIBP :
*--#[ DoPrepareIBP :

**
* A helper procedure that calls PrepareIBP() when previously it has not been
* called.
*
#procedure DoPrepareIBP()
* TODO: reprepare when the input changed.
  #ifndef `IBPYOperator'
    #call PrepareIBP(`prefix'Y)
  #endif
#endprocedure

*--#] DoPrepareIBP :
*--#[ GenerateIBP :

**
* Generates a set of IBP relations. The argument "A" is used for the prefix of
* the (global) expressions holding a relation, Ai == 0. "G" is used for the name
* of integrals. "n" is used for the prefix of index symbols. "D" is the symbol
* for the space-time dimension.
*
* Example:
*   #-
*   Format nospaces;
*   #include ibp.h
*   V k1,k2,p;
*   S D;
*   #define Internals   "k1"
*   #define Externals   "p"
*   #define Propagators "k1.k1,k2.k2"
*   #define Replacements "
*     id k2 = k1 + p;
*   "
*   #call GenerateIBP(A,G,n,D)
*   B G;
*   .sort
*   Collect dum_;
*   P+s;
*   .end
*
* Output:
*   A1=
*      +(-n2)*G(-1+n1,-1+n2)
*      +(-n2)*G(n1,-2+n2)
*      +(-2*n1)*G(-2+n1,n2)
*      +(p.p*n2)*G(n1,-1+n2)
*      +(D)*G(n1,n2)
*      ;
*
*   A2=
*      +(-n2)*G(n1,-2+n2)
*      +(n2-n1)*G(-1+n1,-1+n2)
*      +(-p.p*n2)*G(n1,-1+n2)
*      +(p.p*n1)*G(-1+n1,n2)
*      +(n1)*G(-2+n1,n2)
*      ;
*
* Options:
*   global       - Declare equations as global expressions. (default)
*   local        - Declare equations as local expressions.
*   [no]leeideas - Use only the small set of IBPs. (default: no)
*
* This procedure contains
*   - Declaration statements,
*   - Definition statements,
*   - Executable statements.
*
#procedure GenerateIBP(A,G,n,D,?opts)
  #define prefix "IBPPrivate"

* Internally used functions.
  #define Del   "`prefix'Del"
  #define Del1  "`prefix'Del1"
  #define Delp  "`prefix'Delp"
* Internally used symbols.
  #define xx    "`prefix'xx"
* Wildcards.
  #define x     "`prefix'x"
  #define k1    "`prefix'k1"
  #define k2    "`prefix'k2"
  #define q     "`prefix'q"

* Options.
  #define Decl     "Global"
  #define LeeIdeas "0"
  #do o={`?opts',}
    #ifdef `o'
      #switch `o'
        #case global
          #redefine Decl "Global"
          #break
        #case local
          #redefine Decl "Local"
          #break
        #case leeideas
          #redefine LeeIdeas "1"
          #break
        #case noleeideas
          #redefine LeeIdeas "0"
          #break
        #default
          #message Error: unknown option "`o'"
          #terminate
          #break
      #endswitch
    #endif
  #enddo

  #call DoPrepareIBP()
  #define Y "`IBPYOperator'"

* The next expression number.
  #define N "-1"
  #call GetNextExpr(`A',N)

* Parse Internals.
  #define NI "0"
  #do k={`NormalizedInternals',}
    #ifdef `k'
      #redefine NI "{`NI'+1}"
      #define Internal`NI' "`k'"
    #endif
  #enddo

* Parse Externals.
  #define NE "0"
  #do k={`NormalizedExternals',}
    #ifdef `k'
      #redefine NE "{`NE'+1}"
      #define External`NE' "`k'"
    #endif
  #enddo

* Parse Propagators.
  #define NP "0"
  #do P={`NormalizedPropagators',}
    #ifdef `P'
      #redefine NP "{`NP'+1}"
      #define Propagator`NP' "`P'"
    #endif
  #enddo

  S `n'1,...,`n'`NP',`xx',`x',`x'1,...,`x'`NP';
  V `k1',`k2',`q';
  CF `G',`Del',`Del1',`Delp';

  #define exprs ""

  #if `LeeIdeas' == 0
*   Generate NI*(NI+NE) IBPs.
*     Del(q,k): the derivative of q(mu) with respect to k(mu).
*     Delp(p,k): the scalar product between p(mu) and the derivative of the
*                integral with respect to k(mu).
    #do k={`NormalizedInternals',}
      #ifdef `k'
        #do p={`NormalizedInternals',`NormalizedExternals',}
          #ifdef `p'
            #redefine exprs "`exprs',`A'`N'"
            #redefine NumberOfRelations "{`NumberOfRelations'+1}"
            `Decl' `A'`N++' = `Del'(`p',`k') + `Delp'(`p',`k');
          #endif
        #enddo
      #endif
    #enddo
  #else
*   Generate NI+NE+1 IBPs. See Sect. V of arXiv:0804.3008.
*   (1) NI IBPs.
    #do i=1,`NI'
      #if `i' < `NI'
        #define j "{`i'+1}"
      #else
        #define j "1"
      #endif
      #redefine exprs "`exprs',`A'`N'"
      #redefine NumberOfRelations "{`NumberOfRelations'+1}"
      `Decl' `A'`N++' = `Delp'(`Internal`j'',`Internal`i'');
      #undefine j
    #enddo
*   (2) NE IBPs.
    #do j=1,`NE'
      #redefine exprs "`exprs',`A'`N'"
      #redefine NumberOfRelations "{`NumberOfRelations'+1}"
      `Decl' `A'`N++' = `Delp'(`External`j'',`Internal1');
    #enddo
*   (3) 1 IBP.
    #redefine exprs "`exprs',`A'`N'"
    #redefine NumberOfRelations "{`NumberOfRelations'+1}"
    `Decl' `A'`N++' =
      #do i=1,`NI'
        #define k "`Internal`i''"
        + `Del'(`k',`k') + `Delp'(`k',`k')
        #undefine k
      #enddo
    ;
  #endif

  inexpression `exprs';  * Only for expressons defined in this module.
*   Handle Del(q,k).
    id `Del'(`k1'?,`k1'?) = `D';
    id `Del'(`k1'?,`k2'?) = 0;
*   Convert Delp(p,k) to a sum of
*     Del(p,expr,k): the scalar product between p(mu) and the derivative of expr
*                    with respect to k(mu).
    id `Delp'(`q'?,`k1'?) = 0
      #do i=1,`NP'
        - `n'`i' * `Y'`i'^-1 * `Del'(`q',`Propagator`i'',`k1')
      #enddo
    ;
*   Handle Del(p,expr,k).
    argument `Del',2;
      repeat;
        `NormalizedReplacements';
      endrepeat;
    endargument;
    #do k={`NormalizedInternals',}
      #ifdef `k'
        #do p={`NormalizedInternals',`NormalizedExternals',}
          #ifdef `p'
            id `Del'(`p',`x'?,`k') = `Del1'(`xx'*`x');
            argument `Del1';
              id `xx' * `k'.`k' = 2 * `k'.`p';
              id `xx' * `k'.`q'? = `p'.`q';
              id `xx' = 0;
            endargument;
            id `Del1'(`x'?) = `x';
          #endif
        #enddo
      #endif
    #enddo
*   Multiply G(n1,...,nN) to build the identities.
    multiply `G'(`n'1,...,`n'`NP');
*   Simplifycations.
    #redefine IdentifyYOps "
    id `G'(<`n'1?>,...,<`n'`NP'?>) * <`Y'1^`x'1?>*...*<`Y'`NP'^`x'`NP'?>
      = `G'(<`n'1-`x'1>,...,<`n'`NP'-`x'`NP'>);"
    #call SimplifyIBP()
  endinexpression;
#endprocedure

*--#] GenerateIBP :
*--#[ GetNextExpr :

**
* Gets the number for the next undefined expression.
*
* Example:
*   L Expr1 = 1;
*   L Expr2 = 1;
*   #define N "-1"
*   #call GetNextExpr(Expr,N)
*   L Expr`N' = 1;  * Expr3
*
#procedure GetNextExpr(argA,argDest)
  #do prcI=1,2147483647
    #if exists(`argA'`prcI') == 0
      #redefine `argDest' "`prcI'"
      #breakdo
    #endif
  #enddo
#endprocedure

*--#] GetNextExpr :
*--#[ SimplifyIBP :

**
* An auxiliary procedure for simplifying IBP relations.
* Valid after PrepareIBP() and GenerateIBP() are called.
*
* This procedure contains
*   - Executable statements.
*
#procedure SimplifyIBP()
  `NormalizedReplacements';
  `IdentifySPs';
  repeat;
    `NormalizedReplacements';
  endrepeat;
  `IdentifyYOps';
#endprocedure

*--#] SimplifyIBP :
*--#[ BuildTurbo:
*
#procedure BuildTurbo(TOPO,B,N,file)
*
* Creates the file with the B-equations for use in the reduction algorithm.
* These equations speed up the reductions by a very good factor.
*
* Much of the code here is similar to MakeSimpler.
* Only the eventual id statements need some shifts.
*
#do i = 1,`NUMEQ'
#$done`i' = 0;
#enddo
#create <`file'>
*#write <`file'> "#commentchar %%"
#do j = 1,1
*
* Start with counting:
* The number of terms with leading complexity in $nn
* The number of terms with subleading complexity in $ns
* The total number of terms of the expression in $nt
*
#do i = 1,`NUMEQ'
#$nn`i' = 0;
#$ns`i' = 0;
#$nt`i' = 0;
#if ( `$done`i'' == 0 )
  #$nt`i' = termsin_(`B'`i');
#endif
#enddo
if ( match(Z(-2,?a)) );
#do i = 1,`NUMEQ'
  if ( ( $done`i' == 0 ) && ( expression(`B'`i') ) ) $nn`i' = $nn`i'+1;
#enddo
elseif ( match(Z(-1,?a)) );
#do i = 1,`NUMEQ'
  if ( ( $done`i' == 0 ) && ( expression(`B'`i') ) ) $ns`i' = $ns`i'+1;
#enddo
endif;
B Z;
.sort
#$nnsum = 1;
#$nssum = 1;
#$ntsum = 1;
#do i = 1,`NUMEQ'
  #$nnsum = $nnsum+`$nn`i'';
  #$nssum = $nssum+`$ns`i'';
  #$ntsum = $ntsum+`$nt`i'';
#enddo
#if ( `$nnsum' > 1 )
    #$imin = 0;
    #$nmin = 10000*`$nnsum'+100*`$nssum'+`$ntsum';
    #do i = 1,`NUMEQ'
    #if ( ( `$nn`i'' > 0 ) && ( `$done`i'' == 0 ) )
    #$n = 10000*`$nn`i''+100*`$ns`i''+`$nt`i'';
    #if ( `$n' < `$nmin' )
      #$imin = `i';
      #$nmin = $n;
    #endif
    #endif
    #enddo
*
* At this point we see i as optimal.
* Now the challenge is to write it away as a proper id statement
* to a file that can be used in the eventual reduction program.
*
  #$exp = `B'`$imin';
  #$donefirstinexp = 0;
  #do n = 1,`N'
    #if ( `$nsign`n'' == -1 )
    #redefine extn`n' "?neg0_"
    #else
    #redefine extn`n' "?pos_"
    #endif
    #$plus`n' = 0;
  #enddo
  #do n = `ORDER'
   #if ( `$donefirstinexp' == 0 );
    #inside $exp
     if ( $donefirstinexp == 0 );
    if ( match(Z(-2,?a,n`n'`nsign`n''1,?b)) );
*
*     We want to work with this term
*
      #do nn = 1,`N'
        if ( match(Z(?a,n`nn'`nsign`nn''1,?b)) );
        if ( $nsign`nn' == -1 );
          redefine extn`nn' "?neg_";
          Multiply replace_(n`nn',n`nn'+1);
          $plus`nn' = 1;
        else;
          redefine extn`nn' "?{>1}";
          Multiply replace_(n`nn',n`nn'-1);
          $plus`nn' = -1;
        endif;
        endif;
      #enddo
      $donefirstinexp = 1;
      id rat(x1?$x1,x2?$x2) = 0;
    endif;
     endif;
    #endinside
  #endif
  #enddo
  #$BADEXP = 0;
  #if ( `$donefirstinexp' );
    #inside $exp
* filter turbos that do not decrease complexity
    if (match(Z(-2,?a))) $BADEXP = 1;
    #endinside
  #endif
  #if ( (`$donefirstinexp') && (`$BADEXP' == 0) )
    #inside $exp
      id Z(x?,?a) = Z(?a);
      Multiply replace_(<n1,n1+$plus1>,...,<n`N',n`N'+$plus`N'>);
      id Z(x1?,...,x`N'?) = M`TOPO'/<i1^(x1-n1)>/.../<i`N'^(x`N'-n`N')>;
      multiply 1/M`TOPO';
*     Expand simple rat-functions.
      splitfirstarg rat;
      id rat(0,?xx,y?) = rat(?xx,y);
      id rat(x?,y?number_) = x/y;
      transform rat, addargs(1,last-1);
    #endinside
    #write <`file'> "id,ifmatch->1,M`TOPO'/i1^n1`extn1'%"
    #do nn = 2,`N'
      #write <`file'> "/i`nn'^n`nn'`extn`nn''%"
    #enddo
    #write <`file'> "="
    #write <`file'> "   -M`TOPO'/<i1^n1>/\.../<i`N'^n`N'>*rat(%3$,%3$)*(%3$);", \
      $x2,$x1,$exp
  #endif
    B Z;
    .sort
    #redefine j "0"
  #$done`$imin' = 1;
#endif
#enddo
*#write <`file'> "#commentchar *"
*
#endprocedure
*
*--#] BuildTurbo : 
*--#[ NewCleanup :
*
#procedure NewCleanup(A,B,n,k)
*
* Gaussian elimination for CleanupBasic.
*
B Z;
.sort
#do j = 1,1
*
* Start with counting:
* The number of terms with leading complexity in $nn
* The number of terms with subleading complexity in $ns
* The total number of terms of the expression in $nt
*
#do i = 1,`NUMEQ'
  #$nna`i' = 0;
  #$nnb`i' = 0;
  #$nnc`i' = 0;
  #$nt`i' = termsin_(`A'`i');
#enddo
if ( match(Z(-2,?a)) );
#do i = 1,`NUMEQ'
  if ( expression(`A'`i') ) $nna`i' = $nna`i'+1;
#enddo
elseif ( match(Z(-1,?a)) );
#do i = 1,`NUMEQ'
  if ( expression(`A'`i') ) $nnb`i' = $nnb`i'+1;
#enddo
elseif ( match(Z(0,?a)) );
#do i = 1,`NUMEQ'
  if ( expression(`A'`i') ) $nnc`i' = $nnc`i'+1;
#enddo
endif;
B Z;
.sort
*#do i = 1,`NUMEQ'
*#write <> " `i': `$nna`i'',`$nnb`i'',`$nt`i''"
*#enddo
*
* Count what is still to go. Add one as a nice starting value for minimization
*
#$nnsum = 1;
#$nssum = 1;
#$ntsum = 1;
#do i = 1,`NUMEQ'
  #$nnsum = $nnsum+`$nna`i'';
  #$nssum = $nssum+`$nnb`i'';
  #$ntsum = $ntsum+`$nt`i'';
#enddo
#if ( {`$nnsum'+`$nssum'} > 2 )
*
* Here we still have work to do
* First we have to find out at what 'level'
* We prefer of course complexity 1 with one term only with only one n+1
* The order:
*   One term at subleading and none at leading  -> $level = 1
*   With leading terms  -> $level = 2
*   only subleading terms, more than one  -> $level = 3
*
    #$imin = 0;
  #$level = 10;
    #$nmin = 100*`$nnsum'+10*`$nssum'+`$ntsum';
    #do i = 1,`NUMEQ'
    #if ( `$nna`i'' > 0 )
    #if ( `$level' > 2 )
      #$level = 2;
      #$imin = `i';
      #$nmin = 100*`$nna`i''+10*`$nnb`i''+`$nt`i'';
    #elseif ( `$level' == 2 )
      #$n = 100*`$nna`i''+10*`$nnb`i''+`$nt`i'';
      #if ( `$n' < `$nmin' )
      #$imin = `i';
      #$nmin = $n;
      #endif
    #endif
    #elseif ( `$nnb`i'' == 1 )
    #if ( `$level' > 1 )
      #$level = 1;
      #$imin = `i';
      #$nmin = 10*`$nnb`i''+`$nt`i'';
    #elseif ( `$level' == 1 )
      #$n = 10*`$nnb`i''+`$nt`i'';
      #if ( `$n' < `$nmin' )
      #$imin = `i';
      #$nmin = $n;
      #endif
    #endif
    #elseif ( `$nnb`i'' > 1 )
    #if ( `$level' > 3 )
      #$level = 1;
      #$imin = `i';
      #$nmin = 10*`$nnb`i''+`$nt`i'';
    #elseif ( `$level' == 1 )
      #$n = 10*`$nnb`i''+`$nt`i'';
      #if ( `$n' < `$nmin' )
      #$imin = `i';
      #$nmin = $n;
      #endif
    #endif
    #endif
    #enddo
*
* At this point we see i as optimal.
*

*    #write "Eliminating `$imin' (`$nna`$imin'',`$nnb`$imin'',`$nnc`$imin'',`$nt`$imin'') `$level'"
  #if ( ( `$nna`$imin'' == 0 ) && ( `$nnb`$imin'' == -10 ) )
    #$exp = `A'`$imin';
    #inside $exp
      if ( match(Z(-1,?a)) );
        #do nn = 1,`n'
        #redefine extn`nn' "?"
        if ( match(Z(-1,?a,1`nsign`nn''n`nn',?b)) );
            redefine extn`nn' "`nsign`nn''1";
        endif;
        #enddo
        id Z(-1,?a) = 1;
        id rat(x1?$x1,x2?$x2) = 0;
      endif;
      id  Z(x?,?a) = Z(c+x+1,?a);
    #endinside
    if ( expression(`A'`$imin') ) discard;
    if ( expression(`B'`$imin') == 0 )
      repeat id Z(c?neg0_
        #do nn = 1,`n'
          ,n`nn'`extn`nn''
        #enddo
          ) = -rat($x2,$x1)*(`$exp');
*#write <> "id  Z(c?neg_%"
        #do nn = 1,`n'
*#write <> ",n`nn'`extn`nn''%"
        #enddo
*#write <> ")"
**#write <> "  = -rat(`$x2',`$x1')*(`$exp');"
  #elseif ( ( `$nna`$imin'' == 0 ) && ( `$nnb`$imin'' > 1 ) )
    if ( expression(`A'`$imin') ) Discard;
  #elseif ( `$nna`$imin'' > 1 )
*
*   Now we should look for an integral with ki+1
*
    #$exp = `A'`$imin';
    #$kterm = -1;
    #$ktermn = 0;
*   #do kk = `N',1,-1
*   #do kk = 1,`N'
*   #do kk = {12,14,13,1,11,5,2,6,3,4,7,10,9,8}
*   #do kk = {12,11,13,14,1,5,9,3,7,10,2,4,8,6}
*   #do kk = {1,5,12,11,13,14,9,3,7,10,2,4,8,6}
*   #do kk = {13,14,12,3,2,1,5,10,7,9,4,11,8,6}
*   #do kk = {3,6,11,10,14,13,8,2,1,4,5,12,9,7}
*   #do kk = {1,3,14,4,11,13,6,5,7,8,10,9,12,2}
*   #do kk = {2,4,14,13,6,1,3,5,12,11,8,7,10,9}
*   #do kk = {1,6,11,13,14,12,7,2,9,3,5,10,4,8}
*   #do kk = {1,6,7,14,13,11,12,3,9,2,5,10,4,8}
*   #do kk = {14,13,12,11,10,7,3,9,2,5,4,8,1,6}
*   #do kk = {4,1,14,13,12,3,9,11,10,7,2,5,8,6}
    #do kk = `ORDER'
     #if ( ( `$nsign`kk'' == -1 ) && ( `$ktermn' == 0 ) )
      #inside $exp
      if ( match(Z(-2,?a)) );
        if ( $ktermn == 0 );
          if ( match(Z(-2,?a,n`kk'-1,?b)) );
          $ktermn = 1;
          id  rat(x1?$x1,x2?$x2) = 1;
          $kterm = term_;
          Discard;
          endif;
        endif;
      endif;
      #endinside
     #endif
    #enddo
    #if ( `$ktermn' == 1 )
      if ( expression(`B'`$imin') == 0 );
        id  `$kterm' = -rat(`$x2',`$x1')*(`$exp');
*#write <> "id  `$kterm'"
      endif;
    #else
        #$lhs = FirstBracket_(`A'`$imin');
        #$lhsc = `A'`$imin'[`$lhs'];
        #inside $lhsc;
        id  rat(x1?,x2?) = rat(x2,x1);
        #endinside;
      if ( expression(`B'`$imin') == 0 )
          id  `$lhs' = `$lhs'-`A'`$imin'*`$lhsc';
*#write <> "id  `$lhs'"
    #endif
  #else
      #$lhs = FirstBracket_(`A'`$imin');
      #$lhsc = `A'`$imin'[`$lhs'];
      #inside $lhsc;
      id  rat(x1?,x2?) = rat(x2,x1);
      #endinside;
    if ( expression(`B'`$imin') == 0 )
        id  `$lhs' = `$lhs'-`A'`$imin'*`$lhsc';
*#write <> "   id  `$lhs'"
  #endif
    B Z;
    .sort
    #redefine j "0"
#endif
#enddo
#endprocedure
*
*--#] NewCleanup : 

*--#] proc :
#endif
* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4
