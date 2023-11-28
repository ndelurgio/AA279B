% ------------------------------------------------------------------------------
%
%                           function lambertu
%
%  this function solves the lambert problem for orbit determination and returns
%    the velocity vectors at each of two given position vectors.  the solution
%    uses universal variables for calculation and a bissection technique
%    updating psi.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
% Following line added by Barrows 1/2014
%    mu          - central body gravitational parameter  km^3/sec^2
%    r1          - ijk position vector 1          km
%    r2          - ijk position vector 2          km
%    dm          - direction of motion            'l','s'
%    dtsec       - time between r1 and r2         s
%    nrev        - multiple revoluions            0, 1, ...
%
%  outputs       :
%    v1          - ijk velocity vector            km / s
%    v2          - ijk velocity vector            km / s
%    error       - error flag                     'ok', ...
%
%  locals        :
%    vara        - variable of the iteration,
%                  not the semi-axis
%    y           - area between position vectors
%    upper       - upper bound for z
%    lower       - lower bound for z
%    cosdeltanu  - cosine of true anomaly change  rad
%    f           - f expression
%    g           - g expression
%    gdot        - g dot expression
%    xold        - old universal variable x
%    xoldcubed   - xold cubed
%    zold        - old value of z
%    znew        - new value of z
%    c2new       - c2(z) function
%    c3new       - c3(z) function
%    timenew     - new time                       s
%    small       - tolerance for roundoff errors
%    i, j        - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    dot         - dot product of two vectors
%    findc2c3    - find c2 and c3 functions
%
%  references    :
%    vallado       2001, 459-464, alg 55, ex 7-5
%
% [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec );
% ------------------------------------------------------------------------------

%function [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec );
% Following line commented out by Barrows 1/2014
%function [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec,fid )
% Following line added by Barrows 4/2015 (function renamed; v1_out, v2_out, mu, r1_in, and r2_in added; vo, v, ro, r, and fid removed)
function [v1_out,v2_out,errorl] = AA279lambert_vallado_u ( mu, r1_in,r2_in, dm, nrev, dtsec )

% -------------------------  implementation   -------------------------
% Following line commented out by Barrows 1/2014
%        constmath;
% Following line added by Barrows 1/2014 (twopi added from constmath.m)
        twopi  = 2.0 * pi;
% Following line commented out by Barrows 1/2014 (don't want hardwired mu)
%        constastro;
small = 0.00001; % can affect cases where znew is multiples of 2pi^2
% Following line commented out by Barrows 4/2015 (interplanetary problems were taking 40-60 loops to converge)
%       numiter= 40;
% Following line added by Barrows 4/2015
        numiter= 500;
        errorl  = '      ok';
        psinew = 0.0;

        % Following 2 lines added by Barrows 1/2014
        ro = r1_in;
        r  = r2_in;
        
        magro = mag(ro);
        magr  = mag(r);
        for i= 1 : 3
            vo(i)= 0.0;
            v(i) = 0.0;
        end

        cosdeltanu= nrev + dot(ro,r)/(magro*magr);
        if ( dm == 'l' )  
            vara = -sqrt( magro*magr*(1.0+cosdeltanu) );
        else
            vara =  sqrt( magro*magr*(1.0+cosdeltanu) );
        end
 %fprintf(1,'%11.7f %11.7f nrev %3i %1c \n',cosdeltanu*rad, vara , nrev, dm);

        % ---------------  form initial guesses   ---------------------
        psiold = 0.0;
        psinew = 0.0;
        xold   = 0.0;
        c2new  = 0.5;
        c3new  = 1.0/6.0;

        % --------- set up initial bounds for the bissection ----------
        if ( nrev == 0 )  
            upper=  4.0*pi*pi;
            lower= -4.0*twopi*pi;
            nrev = 0;
        else
            if nrev == 1  
                upper=   16.0*pi*pi; 
                lower=    4.0*pi*pi;   
            else
                upper=  36.0*pi*pi; 
                lower=  16.0*pi*pi;     
            end    
        end

%        chord = sqrt( magro^2 + magr^2 - 2.0*magro*magr*cosdeltanu );
%            nrev = 1;
%        chord = sqrt( magro^2 + magr^2 - 2.0*magro*magr*cosdeltanu );
%        s     = ( magro + magr + chord )*0.5;
%        betam = 2.0* asin( sqrt((s-chord)/chord) );  % comes out imaginary?? jst for hyperbolic??
%        tmin  = ((2.0*nrev+1.0)*pi-betam + sin(betam))/sqrt(mu);

        % -------  determine if  the orbit is possible at all ---------
        if ( abs( vara ) > small )  
            loops  = 0;
            ynegktr= 1;  % y neg ktr
            dtnew = -10.0;
            while ((abs(dtnew-dtsec) >= small) && (loops < numiter) && (ynegktr <= 10))
%       fprintf(1,'%3i  dtnew-dtsec %11.7f yneg %3i \n',loops,dtnew-dtsec,ynegktr );
                if ( abs(c2new) > small )
                    y= magro + magr - ( vara*(1.0-psiold*c3new)/sqrt(c2new) );
                else
                    y= magro + magr;
                end
                % ----------- check for negative values of y ----------
                if (  ( vara > 0.0 ) && ( y < 0.0 ) )  % ( vara > 0.0 ) &
                    ynegktr= 1;
                    while (( y < 0.0 ) && ( ynegktr < 10 ))
                        psinew= 0.8*(1.0/c3new)*( 1.0 ...
                                - (magro+magr)*sqrt(c2new)/vara  );  
                        % -------- find c2 and c3 functions -----------
                        [c2new,c3new] = findc2c3( psinew );
                        psiold = psinew;
                        lower  = psiold;
                        if ( abs(c2new) > small )
                            y= magro + magr - ( vara*(1.0-psiold*c3new)/sqrt(c2new) );
                        else
                            y= magro + magr;
                        end
  %         fprintf(1,'%3i  y %11.7f lower %11.7f c2new %11.7f psinew %11.7f yneg %3i \n',loops,y,lower,c2new,psinew,ynegktr );

                        ynegktr = ynegktr + 1;
                    end % while
                end  % if  y neg

                if ( ynegktr < 10 )  
                    if ( abs(c2new) > small )  
                        xold= sqrt( y/c2new );
                    else
                        xold= 0.0;
                    end
                    xoldcubed= xold*xold*xold;
                    dtnew    = (xoldcubed*c3new + vara*sqrt(y))/sqrt(mu);

                    % --------  readjust upper and lower bounds -------
                    if ( dtnew < dtsec )
                        lower= psiold;
                    end
                    if ( dtnew > dtsec )
                        upper= psiold;
                    end
                    psinew= (upper+lower) * 0.5;

                    % ------------- find c2 and c3 functions ----------
                    [c2new,c3new] = findc2c3( psinew );
                    psiold = psinew;
                    loops = loops + 1;

                    % --- make sure the first guess isn't too close ---
                    if ( (abs(dtnew - dtsec) < small) && (loops == 1) );
                        dtnew= dtsec-1.0;
                    end
                end  % if  ynegktr < 10
%              fprintf(1,'%3i  y %11.7f xold %11.7f dtnew %11.7f psinew %11.7f \n',loops,y,xold,dtnew,psinew );
 %%%             fprintf(1,'%3i  y %11.7f xold %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,xold/sqrt(re),dtnew/tusec,psinew );
%              fprintf(1,'%3i  y %11.7f xold %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,xold/sqrt(re),dtnew/60.0,psinew );
            end % while loop

            if ( (loops >= numiter) || (ynegktr >= 10) )
                errorl= 'gnotconv';
                if ( ynegktr >= 10 )
                    errorl= 'y negati';
                end
            else
                % --- use f and g series to find velocity vectors -----
                f   = 1.0 - y/magro;
                gdot= 1.0 - y/magr;
                g   = 1.0 / (vara*sqrt( y/mu ));  % 1 over g
                for i= 1 : 3
                    vo(i)= ( r(i) - f*ro(i) )*g;
                    v(i) = ( gdot*r(i) - ro(i) )*g;
                end
            end   % if  the answer has converged
        else
% Following line commented out by Barrows 4/2015
%           error= 'impos180';
% Following line added by Barrows 4/2015 ('error' changed to 'errorl')            
            errorl= 'impos180';
        end  % if  var a > 0.0
          
%       fprintf( fid,'psinew %11.5f  %11.5f %11.5f  \n',psinew, dtsec/60.0, xold/rad);
       if errorl ~= '      ok'
 
% Following line commented out by Barrows 1/2014
%           fprintf( fid,'%s ',errorl );
% Following line added by Barrows 1/2014 (fid changed to 1)        
           fprintf( 1,'%s ',errorl ); 
           
       end;
          
% Following 2 lines added by Barrows 1/2014
v1_out = vo';
v2_out = v';

% Following line added by Barrows 4/2015 to allow addition of subfunctions
end % terminates MATLAB function

% Following two subfunctions added by Barrows 4/2015

% ------------------------------------------------------------------------------
%
%                            function mag
%
%  this function finds the magnitude of a vector.  the tolerance is set to
%    0.000001, thus the 1.0e-12 for the squared test of underflows.
%
%  author        : david vallado                  719-573-2600   30 may 2002
%
%  revisions
%    vallado     - fix tolerance to match coe, eq, etc            3 sep 2002
%
%  inputs          description                    range / units
%    vec         - vector
%
%  outputs       :
%    mag         - magnitude
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
% mag = ( vec );
% ----------------------------------------------------------------------------- }

function mag = mag ( vec );

        temp= vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3);

        if abs( temp ) >= 1.0e-16
            mag= sqrt( temp );
          else
            mag= 0.0;
        end
end % terminates MATLAB subfunction

% ------------------------------------------------------------------------------
%
%                           function findc2c3
%
%  this function calculates the c2 and c3 functions for use in the universal
%    variable calculation of z.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    znew        - z variable                     rad2
%
%  outputs       :
%    c2new       - c2 function value
%    c3new       - c3 function value
%
%  locals        :
%    sqrtz       - square root of znew
%
%  coupling      :
%    sinh        - hyperbolic sine
%    cosh        - hyperbolic cosine
%
%  references    :
%    vallado       2001, 70-71, alg 1
%
% [c2new,c3new] = findc2c3 ( znew );
% ------------------------------------------------------------------------------

function [c2new,c3new] = findc2c3 ( znew );

        small =     0.00000001;

        % -------------------------  implementation   -----------------
        if ( znew > small )
            sqrtz = sqrt( znew );
            c2new = (1.0 -cos( sqrtz )) / znew;
            c3new = (sqrtz-sin( sqrtz )) / ( sqrtz^3 );
          else
            if ( znew < -small )
                sqrtz = sqrt( -znew );
                c2new = (1.0 -cosh( sqrtz )) / znew;
                c3new = (sinh( sqrtz ) - sqrtz) / ( sqrtz^3 );
              else
                c2new = 0.5;
                c3new = 1.0 /6.0;
              end
          end
end % terminates MATLAB subfunction

