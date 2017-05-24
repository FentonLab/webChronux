import sys, traceback
import numpy

from numpy import *
from chronuxplus.models import *

import subprocess, datetime, time
import os

from chronuxplus.chronuxObjs import *
from chronuxplus.models import *

import numpy, scipy

import scipy.io
from numpy import *

import subprocess, datetime, time
import os

SECTIONTYPE_TRANSVERSE = "Transverse"
SECTIONTYPE_SAGITTAL = "Sagittal"
SECTIONTYPE_CORONAL = "Coronal"

COLOR_MAP = "jet"

SHOW_SECTIONALS = True
SHOW_SECTIONALS_3D = True

GRID_LIMIT_FRACTION = 0.25
OPACITY_INCREASE = 0.15
OPACITY_DECREASE = 0

COLOR_MAP="jet"
SCALP_RADIUS = 11

NUMBER_SECTIONS = 15
NUMBER_COLUMNS = 5

class ChannelObj(object):
    def __init__(self, ChannelId,  XCoordinate, YCoordinate, ZCoordinate):
        self.channelId  = ChannelId
        self.xCoordinate = XCoordinate
        self.yCoordinate = YCoordinate
        self.zCoordinate = ZCoordinate
    def __unicode__(self):
        return str(self.xCoordinate) + str(self.yCoordinate) + str(self.zCoordinate)

class ParameterObj(object):

    def __init__(self, OpacityIncrease, OpacityDecrease, ColorMap):
        self.show2DSections = ''
        self.show2DTransverseView = ''
        self.show2DCoronalView = ''
        self.show2DSagittalView = ''
        self.show3DSections = ''
        self.numColumns2D = ''
        self.numRows2D = ''
        self.opacityIncrease  = ''
        self.opacityIncrease  = ''
        self.opacityDecrease  = ''
        self.colorMap = ''
    def __unicode__(self):
        return str(self.opacityIncrease) + str(self.opacityDecrease) + str(self.colorMap)

class LbexObj(object):
    def __init__( self):
        self.selectedShell = ''
        self.shells = ''
        self.scalpRadius = ''
        self.scalpConductivity = ''
        self.sphereCenterTolerance = ''
        self.modelRadii = ''
        self.modelConductivity = ''
        self.voxelSize = ''
        self.coordinateOriginX = ''
        self.coordinateOriginY = ''
        self.coordinateOriginZ = ''
        self.seriesCutoff = ''
        self.tailApproxDegree = ''
        self.nMax = ''
        self.concROIVolume = ''
        self.concCutOff = ''
        self.concThreshold = ''
        self.numVoxels = ''
        self.lbexTrialObjs = []
        def __unicode__(self):
            return str ( self.scalpRadius )

def getAnCoefficients(lbexObj):

    # print " in an cn nmax = " + str ( lbexObj.nMax ) 
    
    cn,omegan = getCnCoefficients( lbexObj, lbexObj.nMax)

    # print " in an cn = " + str ( cn ) 

    nk = array ( arange ( lbexObj.seriesCutOff + 1 , lbexObj.nMax ) )

    cO = []
    for index, cnVal in enumerate ( cn ) :
        if omegan[index] != 0:
            cO.append ( cn[index] / omegan[index] )

    cO = cO [ lbexObj.seriesCutOff : lbexObj.nMax-1]

    omegan = omegan [ lbexObj.seriesCutOff : lbexObj.nMax -1]

    rhs=zeros ( lbexObj.tailApproxDegree + 1)
    tmp=zeros(2 * lbexObj.tailApproxDegree + 2);
    nkPower=zeros(len(nk))

    for k in range ( lbexObj.tailApproxDegree + 1 ):
        nkPower = pow ( nk, k )
        tempArray = [x*y for x,y in zip(cO, nkPower) ]
        rhs[k] = sum(  tempArray )
        tempArray = [x/y for x,y in zip(nkPower, omegan) if y!= 0]
        tmp[k] = sum( tempArray )

    tempArray = []

    for k in range ( lbexObj.tailApproxDegree + 2 , 2 * lbexObj.tailApproxDegree + 2):
        tempArray = pow(nk,(k-1)) 

        tempArray = [x/y for x,y in zip(tempArray, omegan) if y!= 0]

        tmp[k-1] = sum ( tempArray ) 

    A = [] 

    for k in range( lbexObj.tailApproxDegree +1) :
        # print " appending = " + str ( tmp[k:k+lbexObj.tailApproxDegree + 1] )  
        A.append ( tmp[k:k+lbexObj.tailApproxDegree + 1] )

    A = numpy.asmatrix(A)

    B = numpy.asmatrix(rhs) 
    B = B.transpose()

    # print " A = " + str ( A ) 
    # print " B = " + str ( B ) 
    an = linalg.inv(A) * B

    return an

def getCnCoefficients(lbexObj, nMax):

    shells = lbexObj.shells

    conductivities = [a.conductivity * float ( lbexObj.scalpConductivity ) for a in shells]
    shellRadii = [ a.radius * float ( lbexObj.scalpRadius ) for a in shells]

    seriesCutoff = lbexObj.seriesCutOff + 1
    conductivity_temp = []              
    conductivity_temp.append(  conductivities[0] ) 
    conductivity_temp.append(  conductivities[2] ) 
    conductivity_temp.append ( conductivities[1] ) 
    conductivity_temp.append ( conductivities[3] ) 

# % These must be normalised wrt scalpRadius
    b = shellRadii[3] / shellRadii[0] 
    b2 = b*b
    c = shellRadii[2] / shellRadii[0]
    c2 = c*c
    d = shellRadii[1] / shellRadii[0]
    d2 = d*d

# conductivity ratios
    conductivityRatio01 = conductivity_temp[0] / conductivity_temp[1] 
    conductivityRatio12 = conductivity_temp[1] / conductivity_temp[2] 
    conductivityRatio23 = conductivity_temp[2] / conductivity_temp[3] 
    
    conductivityRatios = [conductivityRatio01, conductivityRatio12, conductivityRatio23]

    bP = zeros( nMax )
    cP = zeros( nMax )
    dP = zeros( nMax )

    bP[0] = b*b2
    cP[0] = c*c2
    dP[0] = d*d2

    for index in range (1 , nMax-1):

        bP[ index ] = bP[ index - 1 ] * b2
        cP[ index ] = cP[ index - 1 ] * c2
        dP[ index ] = dP[ index - 1 ] * d2

    bP = array(bP)
    cP = array(cP)
    dP = array(dP)

    cn = zeros( nMax )

    for n in range( nMax ):

        np1 = n+2

        x1 = dP[n]* ( bP[n]*(n+1)*(conductivityRatio01-1)*(conductivityRatio12-1)*np1 + cP[n]*(conductivityRatio01*(n+1) + np1)*(conductivityRatio12*(n+1) + np1) )

        x2 = ( (conductivityRatio23*(n+1) + np1) + dP[n]*np1*(conductivityRatio23-1) )

        x3 = np1*cP[n]*( bP[n]*(conductivityRatio01-1)*(conductivityRatio12*np1 + (n+1)) + cP[n]*(conductivityRatio01* (n + 1)  + np1)*(conductivityRatio12-1) )

        x4 = ( (n+1)*(conductivityRatio23-1) + dP[n]*(conductivityRatio23*np1 + (n+1) ) )

        lambda1 = x1 * x2 + x3 * x4

        cn[n] = (  pow (   (2*(n+1)+1 ) , 4  )  )*dP[n]*cP[n]/lambda1

    maxArray = array(arange(1,nMax + 1))
    omegan = [ ( ( b*y ) * (2*y -1)  ) for y in maxArray  ]

    omegan = [x/y for x,y in zip(omegan, bP)]

    return (cn,omegan)

def getLegendre ( xTemp, degree) :

    x2 = [x*x for x in xTemp]
    p0 = [ ( x, (1.5*y - 0.5 ) , x * (2.5 * y -1.5 )  ) for x,y in zip ( xTemp, x2 ) ]
    p1 = [ ( -sqrt ( 1 - y ) , -3* x * sqrt ( 1 - y ) , -sqrt ( 1 - y ) * (2.5 * y -1.5 )  ) for x,y in zip ( xTemp, x2 ) ]

    return p0, p1

def getEEGDipoleInSphere3D(eegSensorLocations, dipoleLocations, lbexObj):

    # an = getAnCoefficients ( NMAX, TAIL_APPROX_DEGREE , SERIES_CUTOFF, SCALP_RADIUS, conductivities, shellRadii) 
    an = getAnCoefficients ( lbexObj) 

    an = array ( an.tolist() ) 

    # [ cn, junk ] = getCnCoefficients( SERIES_CUTOFF + 1 , SCALP_RADIUS, conductivities , shellRadii)

    nMax = lbexObj.seriesCutOff + 1

    [ cn, junk ] = getCnCoefficients( lbexObj, nMax)

    cn = cn[:len(cn)-1]

    ns = array ( arange ( 1,lbexObj.seriesCutOff +1) )

    # print " lbexObj.seriesCutOff = " + str ( lbexObj.seriesCutOff ) 
    # print " cn = " + str (  cn ) 
    # print " shape cn = " + str ( shape ( cn ) ) 
    # print " an = " + str ( an ) 
    # print " ns = " + str ( ns ) 

    cnt = cn - ( an[0] + ns*( an[1] + ns*( an[2] + ns*an[3] ) ) )

    const = 1 / (  4 * math.pi * lbexObj.scalpConductivity * lbexObj.scalpRadius * lbexObj.scalpRadius )

    cnt = const * cnt; 
    cn = const*cn; 
    an = const*an

    tolerance = lbexObj.sphereCenterTolerance * lbexObj.scalpRadius

    rd2 = array ( [x[0]*x[0] + x[1]*x[1] + x[2]*x[2] for x in dipoleLocations] ) 

    inc = array ( [ index for index, x in enumerate(rd2) if x > tolerance*tolerance ] )  
    
    dIndx = array ( [[1+3*inc ],  [2+3*inc], [3+3*inc] ] ) 
    dIndx = dIndx.flatten(1)

    nwDipoleLocations = []
    nwrd2 = []

    for index in inc:
        nwDipoleLocations.append( dipoleLocations [ index] )
        nwrd2.append ( rd2[ index])

    nSensors = len ( eegSensorLocations )
    nDipoles = len ( nwDipoleLocations)

    rd = [ math.sqrt( a ) for a in nwrd2 ]
    r0 = [ (x[0]/y, x[1]/y, x[2]/y) for x,y in zip ( nwDipoleLocations, rd) if y!= 0 ] 

    f = [x/lbexObj.scalpRadius for x in rd]

    s0 = [  ( a.xCoordinate / math.sqrt (a.xCoordinate * a.xCoordinate + a.yCoordinate*a.yCoordinate + a.zCoordinate * a.zCoordinate ) , a.yCoordinate / math.sqrt (a.xCoordinate * a.xCoordinate + a.yCoordinate*a.yCoordinate + a.zCoordinate * a.zCoordinate ) , a.zCoordinate / math.sqrt (a.xCoordinate * a.xCoordinate + a.yCoordinate*a.yCoordinate + a.zCoordinate * a.zCoordinate ) ) for a in eegSensorLocations]
    
    leadField = []

    for ns in range ( nSensors ) :

        temp =  tile ( eegSensorLocations[ns], ( nDipoles, 1 ) )

        rds = [ ( x[0] * y[0].xCoordinate + x[1] * y[0].yCoordinate + x[2] * y[0].zCoordinate ) for x,y in zip ( nwDipoleLocations,  temp ) ]
        
        t01 = [ ( x[0] * y, x[1] * y, x[2] * y ) for x,y in zip ( nwDipoleLocations, rds  ) ]

        t02 = [ ( x[0].xCoordinate * y, x[0].yCoordinate * y, x[0].zCoordinate * y ) for x,y in zip ( temp, nwrd2 )  ] 

        t0 = [ ( x[0]-y[0], x[1]-y[1], x[2]-y[2] ) for x,y in zip ( t02, t01 )  ] 

        tmp = [ math.sqrt ( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ) for x in t0 ] 

        iZeros = [ index for index, x in enumerate ( tmp ) if x > tolerance ]

        for index in iZeros:
            t0[ index ] = ( t0 [ index ][0]   / tmp [ index ] , t0 [ index ][1]   / tmp [ index ], t0 [ index ][2]   / tmp [ index ] )  # Unit tangent vector # nDipoles x 3

        s0Temp = tile ( s0[ns] , ( nDipoles,1) ) 

        xTemp  = [ ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] ) for x,y in zip ( r0,s0Temp)  ]
        xTemp2  = [ pow( ( x ),2) for x in xTemp ]

        xmf = [ ( x-y ) for x,y in zip ( xTemp, f ) ]

        r2i = [ 1/ ( 1 - x * ( y + z)  ) for x,y,z in zip ( f, xTemp, xmf) ]

        ri = sqrt( r2i)

        tmp = an[3]*r2i

        aL = [ ( x + y* ( ( 5*z-4) + y* ( x * ( z - 9 ) + y * (  ( 10 - z - z ) - y * ( x+y) ) ) ) )*k  for x,y,z,k in zip ( xTemp, f, xTemp2, tmp ) ]

        aL = [ ( (an[2]*( x*( 1 + x*y ) - y*( 2 - y*z ) ) + k )* m ) [0]   for x,y,z,k,m in zip ( xTemp, f, xmf, aL,r2i ) ]

        aL = [ ((  an[1]*x + y )* z ) [0]  for x,y,z in zip ( xmf, aL,r2i ) ]

        aL = [ ((  (an[0]*( 1 - 1/x ) / y) + z ) * x)[0]  for x,y,z in zip ( ri, f, aL) ]

        aM = [ ( 1 + y*( 5*x + y*( (z-10) - y*( k - y - y - y ) ) ) ) * m for x,y,z,k,m in zip ( xTemp,f,xTemp2,xmf,tmp) ]

        aM = [ ( ( an[2]*( 1 + y*( x - y ) ) + z) * k ) [0] for x,y,z,k in zip ( xmf,f,aM,r2i) ]

        aM = [ ((  an[1] + x ) * y ) [0] for x,y in zip ( aM, r2i) ]

        aM = [ ( ( (an[0]*(x + 1) / ( x - y*z*x + 1 )) + k ) *x* sqrt( 1 - m )  ) [0] for x,y,z,k,m in zip ( ri, f, xTemp, aM, xTemp2) ]

        [ p0, p1 ] = getLegendre( xTemp, lbexObj.seriesCutOff ) 

        p1 = [ (x[0] / 1, x[1] / 2, x[2] / 3 ) for x in p1 ] 

        tmpLf =  [ ( cnt[0])*x[0] + y*(     cnt[1]*x[1] +           cnt[2]*y*x[2] ) + z for x,y,z in zip ( p0, f, aL ) ]

        tmpLf2 = [ cnt[0] * x[0] + y*( 0.5*cnt[1]*x[1] + 0.3333333*cnt[2]*y*x[2] ) + z for x,y,z in zip ( p1, f,  aM) ]

        tmpLfFinal = array ( [ ( k[0]* x + z[0]*y , k[1]* x + z[1]*y , k[2]* x + z[2]*y  )   for x,y,z, k  in zip ( tmpLf, tmpLf2, t0, r0 ) ] )

        tmpLfFinal = reshape(tmpLfFinal, 3*nDipoles)

        leadField.append ( tmpLfFinal ) 

    lfTranspose = array(leadField)

    if  nDipoles != len( dipoleLocations): 
        dIndx = []
        inc = [ index for a, index in enumerate ( rd2 ) if a <= tolerance*tolerance ]
        inc=inc.transpose()
        for i in range(1,3):
            dIndx.append( array[i+3*(inc-1) ] ) 
            dIndx = dIndx.flatten(1) 
        leadfield[dIndx] = cn[0]*eegSensorLocations / lbexObj.scalpRadius

    return lfTranspose

def createleadField ( lbexObj, sphereElectrodeLeadList ):

    #sphereElectrodeLeadList = leadList

    #lbexAnalysisDetail = LBEXAnalysisDetail.objects.filter ( analysisDetail = analysisDetail ) [0]

    #sphereElectrodeLeadList = analysisDetail.leadList
    
    print ("************** " + str(lbexObj.refLeadListChannel))

    sphereElectrodeHeadBox = sphereElectrodeLeadList.headbox

    sphereHeadBoxChannels = []

    leadListChannelMap = {}

    leadListChannels = LeadListChannel.objects.filter ( leadList = sphereElectrodeLeadList ) 

    useEGI1020Mapping = lbexObj.useEGI1020Mapping
    refLeadListChannel = lbexObj.refLeadListChannel

    egiLeadListChannelList = []
    ten20ChannelLeadList = []

    if useEGI1020Mapping:
        ten20Channels = EGI1020Assignment.objects.all()
        for ten20Channel in ten20Channels:
            egiLeadListChannelList.append ( ten20Channel.egiLeadListChannel ) 
            ten20ChannelLeadList.append ( ten20Channel.tenTwentyChannel ) 

    for leadListChannel in leadListChannels:

        if useEGI1020Mapping:
            if leadListChannel not in egiLeadListChannelList:
                continue
        # print " useEGI1020Mapping = " + str ( useEGI1020Mapping ) + " channel = " + str ( leadListChannel.headBoxChannel ) 
        #print " leadListChannel = " + str ( leadListChannel.id ) 
        EGI1020AssignmentChannel = EGI1020Assignment.objects.filter ( egiLeadListChannel = leadListChannel )[0]
        tenTwentyChannel = EGI1020AssignmentChannel.tenTwentyChannel
        headBoxChannel = leadListChannel.headBoxChannel
        EGI1020AssignmentChannel = EGI1020Assignment.objects.filter ( egiLeadListChannel = leadListChannel )[0]
        tenTwentyChannel = EGI1020AssignmentChannel.tenTwentyChannel
        headBoxChannel = tenTwentyChannel.headBoxChannel
        # headBoxChannel = leadListChannel.headBoxChannel
        if headBoxChannel.xCoordinate == 0 and headBoxChannel.yCoordinate == 0 and headBoxChannel.zCoordinate == 0:
            continue
        sphereHeadBoxChannels.append(headBoxChannel)
        
    channelDistances = array ( [ sqrt(a.xCoordinate * a.xCoordinate + a.yCoordinate * a.yCoordinate + a.zCoordinate * a.zCoordinate) for a in sphereHeadBoxChannels] ) 
    sensorLocations =  [ChannelObj( a.channelId, float (lbexObj.scalpRadius) * ( a.xCoordinate/float(b) ), float (lbexObj.scalpRadius) * ( a.yCoordinate/float(b) ) , float (lbexObj.scalpRadius) * ( a.zCoordinate/float(b) )  )  for a,b in zip ( sphereHeadBoxChannels, channelDistances) if float ( b ) != 0 ]

    #print " sensorLocations " + str ( sensorLocations ) 
    refHeadboxChannels = []

    # refHeadBoxChannel = HeadBoxChannel.objects.filter ( headbox = sphereElectrodeHeadBox, channelId = REF_HEADBOX_CHANNEL) [0]
    refLeadListEGIChannel = EGI1020Assignment.objects.filter ( egiLeadListChannel = lbexObj.refLeadListChannel ) [0]
    refHeadBoxChannel = refLeadListEGIChannel.tenTwentyChannel.headBoxChannel
    refHeadboxChannels.append ( refHeadBoxChannel ) 

    channelDistances = array ( [ sqrt(a.xCoordinate * a.xCoordinate + a.yCoordinate * a.yCoordinate + a.zCoordinate * a.zCoordinate) for a in sphereHeadBoxChannels] ) 
    sensorLocations =  [ChannelObj( a.channelId, float (lbexObj.scalpRadius) * ( a.xCoordinate/float(b) ), float (lbexObj.scalpRadius) * ( a.yCoordinate/float(b) ) , float (lbexObj.scalpRadius) * ( a.zCoordinate/float(b) )  )  for a,b in zip ( sphereHeadBoxChannels, channelDistances) if float ( b ) != 0 ]

    refChannelDistances = array ( [ sqrt(a.xCoordinate * a.xCoordinate + a.yCoordinate * a.yCoordinate + a.zCoordinate * a.zCoordinate) for a in refHeadboxChannels] ) 
    refSensorLocations =  [ChannelObj( a.channelId, float (lbexObj.scalpRadius) * ( a.xCoordinate/float(b) ), float (lbexObj.scalpRadius) * ( a.yCoordinate/float(b) ) , float (lbexObj.scalpRadius) * ( a.zCoordinate/float(b) )  )  for a,b in zip ( refHeadboxChannels, channelDistances) if float ( b ) != 0]

    # shells = array(Shell.objects.all())

    # print " sensorLocations " + str ( sensorLocations ) 

    xCoords = [ a.xCoordinate for a in sensorLocations ] 
    yCoords = [ a.yCoordinate for a in sensorLocations ] 
    zCoords = [ a.zCoordinate for a in sensorLocations ] 

    xLimits =  [floor( min ( xCoords ) ) , ceil ( max ( xCoords )  ) ] 
    yLimits =  [floor( min ( yCoords ) ) , ceil ( max ( yCoords )  )]  
    zLimits =  [floor( min ( zCoords ) ) , ceil ( max ( zCoords )  )]  
  
    gridLimit = [ xLimits, yLimits, zLimits ]

    fraction = 0.25 # between 0 & 0.5

    fvs = fraction * lbexObj.voxelSize

    xGrid = arange ( xLimits[0]+fvs ,   xLimits[1]-fvs, lbexObj.voxelSize )
    yGrid = arange ( yLimits[0]+fvs ,   yLimits[1]-fvs, lbexObj.voxelSize )
    zGrid = arange ( zLimits[0]+fvs ,   zLimits[1]-fvs, lbexObj.voxelSize )

    nCubes = len(xGrid) * len(yGrid) * len(zGrid)

    x = zeros(nCubes)
    y = zeros(nCubes)
    z = zeros(nCubes)
    
    inside = zeros(nCubes)

    xArray = []
    yArray = []
    zArray = []

    [x,y,z] = numpy.meshgrid( xGrid  , yGrid , zGrid )
    [x,y,z] = meshgrid( xGrid  , yGrid , zGrid )
    xArray = x.flatten(1)
    yArray = y.flatten(1)
    zArray = z.flatten(1)

    radialDistanceArray = [ a*a + b*b + c*c for a,b,c in zip ( xArray, yArray, zArray ) ]

    inside = [index for index,elem in enumerate(radialDistanceArray) if elem<=lbexObj.scalpRadius * lbexObj.scalpRadius]

    grid = [ (x,y,z) for x, y, z in zip ( xArray, yArray, zArray ) ]

    dipoleLocations = []

    for index in inside:
        dipoleLocations.append( grid [ int(index)] )

    numVoxels = len(inside)

    eegLeadField = getEEGDipoleInSphere3D(sensorLocations, dipoleLocations, lbexObj)

    referencedEEGLeadField = numpy.asmatrix(eegLeadField)
            
    return sensorLocations, numVoxels, referencedEEGLeadField , dipoleLocations

#def voxelNeighbours( roiVolume, roiCenteredOnVoxel, voxelCentroids, numVoxels ):

    #refEegLeadField = getEEGDipoleInSphere3D( refSensorLocations, dipoleLocations, lbexObj )

    #referencedEEGLeadField = []

    #for sensorIndex in range ( len ( sensorLocations ) ) :
        #referencedEEGLeadFieldElem  = [ ( x - y ) for x,y in zip ( eegLeadField [ sensorIndex ] , refEegLeadField[0] ) ]
        #referencedEEGLeadField.append ( referencedEEGLeadFieldElem ) 

    #referencedEEGLeadField = numpy.asmatrix(referencedEEGLeadField)
            
    #return numVoxels, referencedEEGLeadField , dipoleLocations

def voxelNeighbours( roiVolume, roiCenteredOnVoxel, voxelCentroids, numVoxels ):

    # print " numVoxels = " + str ( numVoxels ) + "  roiCenteredOnVoxel = " + str ( roiCenteredOnVoxel ) + " voxelCentroids = " + str ( voxelCentroids [roiCenteredOnVoxel] ) + " shape voxelCentroids = " + str ( shape ( voxelCentroids ) ) 

    roiColumns = '' 
    roiVoxelIndices = '' 

    centerVoxel = voxelCentroids[ roiCenteredOnVoxel]

    #print " numVoxels = " + str ( numVoxels ) + " roiVolume = " + str ( roiVolume )
    #print " type numVoxels = " + str ( type ( numVoxels ) ) + " type  roiVolume = " + str ( type ( roiVolume ) )

    #print " ceil = " + str ( ceil(float ( roiVolume ) * float ( numVoxels ) ) ) 

    numROI = min( numVoxels, ceil(float ( roiVolume ) * float ( numVoxels )) )

    centerVoxel = voxelCentroids[ roiCenteredOnVoxel]

    numROI = min( numVoxels, ceil(roiVolume * numVoxels) )
    centerVoxelTile = tile( centerVoxel, ( numVoxels, 1 )) 
    diffVoxels = [ (  x[0] - y[0] , x[1] - y[1] , x[2] - y[2] ) for x,y in zip ( voxelCentroids, centerVoxelTile ) ]
    sumDiffVoxels = [ (  x[0] * x[0] + x[1] * x[1] + x[2] * x[2] ) for x in diffVoxels ]

    argSortVoxels = numpy.argsort ( sumDiffVoxels ) 

    roiVoxelIndices = argSortVoxels [ : numROI ] 

    roiColumns = array ( [ ( 3*x , 1+3*x  , 2+3*x  )  for x in roiVoxelIndices ] ) 

    roiColumns = roiColumns.flatten()

    return roiColumns, roiVoxelIndices

# def ddConcentration( kernel,  voxelCentroids, lbexObj, lbexAnalysisDetail ):
def ddConcentration( eegLeadField,  dipoleLocations, lbexObj, numTapers, concEigenVectorFilePath, concEigenValueFilePath, saveVoxelInfo ):

    try: 

        # lbexFilePath = lbexAnalysisDetail.analysisDetail.lbexFilePath
        # numTapers = lbexAnalysisDetail.analysisDetail.numTapers

        voxelInfoObjList = []

        numVoxels = len( dipoleLocations )
        numChannels = len( eegLeadField)

        szMin = min ( len( eegLeadField) , len( eegLeadField[0]) ) 

        dipoleLocations0 =  dipoleLocations[0] 
        dipoleLocationsTile =  tile ( dipoleLocations[0] , ( len ( dipoleLocations ) -1, 1) ) 

        dipoleLocationsRest =  dipoleLocations[1:] 

        diffVoxels = [ (  x[0] - y[0] , x[1] - y[1] , x[2] - y[2] ) for x,y in zip ( dipoleLocationsTile , dipoleLocationsRest ) ]

        sumDiffVoxels = [ x[0] * x[0] + x[1] * x[1] + x[2] * x[2] for x in diffVoxels ]

        minSumDiffVoxels = min ( sumDiffVoxels ) 

        voxelLengthScale = sqrt ( minSumDiffVoxels ) 

        sourceSpaceVolume = numVoxels * pow ( ( voxelLengthScale ) , 3 ) 

        V, S, U = linalg.svd( eegLeadField, full_matrices=False )
        V, S, U = linalg.svd( kernel, full_matrices=False )

        totalS = sum(S)

        rsum = 0
        i = 0

        #lbexFilePath = lbexAnalysisDetail.analysisDetail.lbexFilePath

        while rsum < lbexObj.concCutOff:
            rsum = rsum + S[i]/totalS 
            i = i+1

        cutIndx = i
        S = S[:i]
        sInv = [1/x for x in S]
        diagSInv = diag(sInv)

        cutUS = numpy.matrix ( U[ : cutIndx ].transpose() ) * numpy.matrix ( diagSInv ) 

        # concEigenVectorFilePath = lbexAnalysisDetail.concEigenVectorFilePath
        # concEigenValueFilePath = lbexAnalysisDetail.concEigenValueFilePath

        #print " concEigenVectorFilePath = " + str ( concEigenVectorFilePath )
        #print " concEigenValueFilePath = " + str ( concEigenValueFilePath )
        # concEigenValueFilePath = lbexAnalysisDetail.concEigenValueFilePath

        if concEigenVectorFilePath != '':
            cvFile = file ( concEigenVectorFilePath , "a")
        if concEigenValueFilePath != '':
            ceFile = file( concEigenValueFilePath , "a")

        for voxelNum in range ( numVoxels , 0 , -1 ) :

            if voxelNum == 1: 
                break
            [ roiColumns, roiVoxelIndices ]  = voxelNeighbours( lbexObj.concROIVolume, voxelNum - 1, dipoleLocations, numVoxels )
        concEigenVectorFilePath = lbexAnalysisDetail.concEigenVectorFilePath
        concEigenValueFilePath = lbexAnalysisDetail.concEigenValueFilePath

        cvFile = file ( concEigenVectorFilePath , "a")
        ceFile = file( concEigenValueFilePath , "a")

        for voxelNum in range ( numVoxels , 0 , -1 ) :

            [ roiColumns, roiVoxelIndices ]  = voxelNeighbours( lbexObj.concROIVolume, voxelNum - 1, voxelCentroids, numVoxels )

            svdV = V[roiColumns]

            [ Up, E, Vp ] =linalg. svd( svdV , full_matrices=False)

            Vp = Vp.transpose()

            powE = pow ( E , 2) 
            total = sum(powE)

            rsum = 0
            mStar = 0
            while rsum < lbexObj.concThreshold:
                rsum = rsum + powE[ mStar ]/total
                mStar = mStar +1

            concVectors = numpy.matrix ( cutUS ) * numpy.matrix ( Vp[ :cutIndx, :mStar ] ) 

            cV = concVectors.transpose()
            cE = tile( powE[ :mStar ], (1, numTapers ) )
            # cEDisplay = tile( powE[ :mStar ], (1, 1 ) )
            cEFirst = cE[0]

            #print " CE size = " + str ( shape ( cE ) ) 
            #print " CE first size = " + str ( shape ( cEFirst [:mStar] ) ) 
            
            # print " for voxel = " + str ( voxelNum )  + " cutUS = " + str ( cutUS ) + " cutIndx = " + str ( cutIndx ) + " mStar = " + str ( mStar )
            if saveVoxelInfo :
                voxelInfo = VoxelInfo ( lbexAnalysisDetail = lbexAnalysisDetail, voxelNum = voxelNum, numEigenValues = mStar) 
                voxelInfo.save()
            # print " for voxel = " + str ( voxelNum )  + " cutUS = " + str ( cutUS ) + " cutIndx = " + str ( cutIndx ) + " mStar = " + str ( mStar )
            voxelInfo = VoxelInfo ( lbexAnalysisDetail = lbexAnalysisDetail, voxelNum = voxelNum, numEigenValues = mStar) 
            voxelInfo.save()

            # print " for voxel = " + str ( voxelNum )  + " cutUS = " + str ( cutUS ) + " cutIndx = " + str ( cutIndx ) + " mStar = " + str ( mStar )
            # print " cV = " + str ( cV ) 
            # print " cE = " + str ( tile( E[ :mStar ], (1, 1) ) ) 
            if concEigenVectorFilePath != '':
                numpy.savetxt ( cvFile, cV, delimiter = ",")
            if concEigenValueFilePath != '':
                numpy.savetxt ( ceFile, tile( powE[ :mStar ], (1, 1) ), delimiter = ",")

            voxelInfoObj = VoxelInfoObj ( )

            voxelInfoObj.voxelNum = voxelNum
            voxelInfoObj.numEigenVectors = len ( cV )
            voxelInfoObj.concEigenVectors = cV
            voxelInfoObj.concEigenValues = cEFirst [:mStar]

            voxelInfoObjList.append ( voxelInfoObj ) 

        if concEigenVectorFilePath != '':
            cvFile.close()
        if concEigenValueFilePath != '':
            ceFile.close()

            #return voxelInfoObjList

            numpy.savetxt ( cvFile, cV, delimiter = ",")
            numpy.savetxt ( ceFile, tile( powE[ :mStar ], (1, 1) ), delimiter = ",")

        cvFile.close()
        ceFile.close()

    except:

        traceback.print_exc(file=sys.stdout) 

def computeLbexAllData( analysisDetail, submittedJob, machines, startJobStatusCode, endJobStatusCode, datafile, jobProcessLbexType ):
    
      lbexAnalysisDetail = LBEXAnalysisDetail.objects.filter ( analysisDetail = analysisDetail ) [0]
      datafile = analysisDetail.datafile

      machine = ""
      for m in machines:
          if m.isMaster:
              machine = m
              break

      jobProcessList = []

      commandString = machine.mappedBasePath + "workspace/Chronux++/lbex"
      command = [ commandString ]
      arguments = []

      arguments.append(str(datafile.id))
      arguments.append(str(analysisDetail.id))
      arguments.append(str(machine.name))
      arguments.append(str(machine.metadataDatabase.name))
      
      dbUrl = str(machine.metadataDatabase.ipAddress) + ":" + str(machine.metadataDatabase.port)
      arguments.append(dbUrl)
      arguments.append(str(machine.metadataDatabase.user))
      arguments.append(str(machine.metadataDatabase.password))
      arguments.append(str(datafile.datafileType))
      # arguments.append(str(lbexPointEvent.id))

      command.extend(arguments)

      #print " analysis command = " + str(command)
      ps = subprocess.Popen( command , stdin = subprocess.PIPE, stdout = subprocess.PIPE )
      #print " ps id for machine " + machine.name + " is " + str(ps) 

      jobProcess = JobProcess(name = "Processing for source localisation of " + str(datafile), jobProcessType = jobProcessLbexType, machine = machine, processId = ps.pid, jobStatusCode = startJobStatusCode, analysisDetail = analysisDetail, submittedJob = submittedJob, startedTime = datetime.datetime.now()) 
      jobProcess.save()
      jobProcessList.append ( jobProcess )

      ps.communicate()

      jobProcess.jobStatusCode = endJobStatusCode
      jobProcess.completedTime = datetime.datetime.now()
      jobProcess.save()    
 
def updateLbexPointEvents (lbexAnalysisDetail) :

    # lbex point events 

      mongoClient = MongoClient()
      pointEventMongoDB = mongoClient.chronuxMongoDB
      pointEventMongoDBFile = pointEventMongoDB.pointEventMongoDBFile
      lbexPointEventMongoDBFile = pointEventMongoDB.lbexPointEventMongoDBFile

      lbexPointEvents = LBEXPointEvent.objects.filter ( lbexAnalysisDetail = lbexAnalysisDetail ) 

      for lbexPointEvent in lbexPointEvents:

          # print "for lbex point event = " + str ( lbexPointEvent.id ) 

          lbexPointEventDetails = lbexPointEventMongoDBFile.find( { "lbexPointEventId" : lbexPointEvent.id } )

          for lbexPointEventDetail in lbexPointEventDetails:

            # print "for lbex point event detail = " + str ( lbexPointEventDetail["pointEventIdentifier"] ) 

            pointEventIdentifier = lbexPointEventDetail["pointEventIdentifier"]
            data = pointEventIdentifier.split("-")
            pointEventId = data[0]
            startPoint = data[1]

            pointEvent = PointEvent.objects.get ( pk = int(pointEventId) )

            pointEventDetailList = pointEventMongoDBFile.find( { "pointEventId" : int(pointEventId), "startPoint": str(startPoint) } ) 

            # print "for lbex point event detail = " + str ( lbexPointEventDetail["pointEventIdentifier"] ) 

            for pointEventDetail in enumerate ( pointEventDetailList ) : 
            # print "startPoint = " + pointEventDetail["startPoint"] + " epochDescription = " + pointEventDetail["epochDescription"] + " eventNum = " + pointEventDetail["eventNum"]
                pointEventDetailObjs = PointEventDetail.objects.filter ( pointEvent = pointEvent, startPoint = startPoint ) 

                if len ( pointEventDetailObjs ) > 0:
                    pointEventDetail = pointEventDetailObjs[0]
                    lbexPointEventDetail = LBEXPointEventDetail(lbexPointEvent = lbexPointEvent, pointEventDetail = pointEventDetail )

                    # print "adding lbex point event detail = " + str ( startPoint ) 

                    lbexPointEventDetail.save()

          lbexPointEventDetailRecs =  LBEXPointEventDetail.objects.filter (lbexPointEvent = lbexPointEvent)

          lbexPointEvent.numLbexPointEvents = len ( lbexPointEventDetailRecs )
          lbexPointEvent.save()

def computeLbexPointEvent( analysisDetail, submittedJob, machines, startJobStatusCode, endJobStatusCode, datafile, jobProcessLbexType ):
    
      lbexAnalysisDetail = LBEXAnalysisDetail.objects.filter ( analysisDetail = analysisDetail ) [0]

      datafile = analysisDetail.datafile

      jobProcessList = []

      machine = ""
      for m in machines:
          if m.isMaster:
              machine = m
              break

      lbexPointEvents  = LBEXPointEvent.objects.filter ( lbexAnalysisDetail = lbexAnalysisDetail ) 
    
      for lbexPointEvent in lbexPointEvents:

        commandString = machine.mappedBasePath + "workspace/Chronux++/lbexp"
        command = [ commandString ]
        arguments = []

        arguments.append(str(datafile.id))
        arguments.append(str(analysisDetail.id))
        arguments.append(str(machine.name))
        arguments.append(str(machine.metadataDatabase.name))

        dbUrl = str(machine.metadataDatabase.ipAddress) + ":" + str(machine.metadataDatabase.port)
        arguments.append(dbUrl)
        arguments.append(str(machine.metadataDatabase.user))
        arguments.append(str(machine.metadataDatabase.password))
        arguments.append(str(datafile.datafileType))
        arguments.append(str(lbexPointEvent.id))

        command.extend(arguments)

        #print " analysis command = " + str(command)
        ps = subprocess.Popen( command , stdin = subprocess.PIPE, stdout = subprocess.PIPE )
        #print " ps id for machine " + machine.name + " is " + str(ps) 

        jobProcess = JobProcess(name = "Processing lbex point event " + str ( lbexPointEvent.id ) + "  for source localisation of " + str(datafile), jobProcessType = jobProcessLbexType, machine = machine, processId = ps.pid, jobStatusCode = startJobStatusCode, analysisDetail = analysisDetail, lbexPointEvent = lbexPointEvent, submittedJob = submittedJob, startedTime = datetime.datetime.now()) 
        jobProcess.save()
        jobProcessList.append ( jobProcess )

        ps.communicate()

        jobProcess.jobStatusCode = endJobStatusCode
        jobProcess.completedTime = datetime.datetime.now()
        jobProcess.save()    
 
   # insert records of lbex point event and point event details into mysql DB for Mysql query for display
 
      jobProcess = JobProcess(name = "Adding database records for point event " + str ( lbexPointEvent.id ) + "  for source localisation of " + str(datafile), jobProcessType = jobProcessLbexType, machine = machine, processId = ps.pid, jobStatusCode = startJobStatusCode, analysisDetail = analysisDetail, lbexPointEvent = lbexPointEvent, submittedJob = submittedJob, startedTime = datetime.datetime.now()) 
      jobProcess.save()
      jobProcessList.append ( jobProcess )
                
    # lbexTrialObjs = lbexObj.lbexTrialObjs

    # for lbexTrialObj in lbexTrialObjs:

    #     lbexPointEventId = lbexTrialObj.lbexPointEventId
    #     lbexPointEvent = LBEXPointEvent.objects.get ( pk = lbexPointEventId ) 

    #     pointEventDetailObjs = lbexTrialObj.pointEventDetailObjs

    #     for pointEventDetailObj in pointEventDetailObjs:
    #         lbexPointEventDetail = LBEXPointEventDetail(lbexPointEvent = lbexPointEvent, startPoint = int ( pointEventDetailObj["startPoint"]), epochDescription = pointEventDetailObj["epochDescription"], eventNum = int ( pointEventDetailObj["eventNum"] ) )
    #         lbexPointEventDetail.save()

      jobProcess.jobStatusCode = endJobStatusCode
      jobProcess.completedTime = datetime.datetime.now()
      jobProcess.save()    

def computeSourceLocalization ( analysisDetail, submittedJob, machines, startJobStatusCode, endJobStatusCode, datafile, jobProcessLbexType, epochType, performLbexAnalysis ):
# def computeSourceLocalization ( ):
    # analysisDetail = AnalysisDetail.objects.get ( pk = 8 ) 

    pslist = []
    jobProcessList = []
    numJobsCompleted = 0

    subtractFromBaseLine = analysisDetail.subtractFromBaseLine
    
    submittedJob.jobStatusCode = startJobStatusCode
    submittedJob.save()

    masterMachine = ""
    for machine in machines:
        if machine.isMaster:
                masterMachine = machine
                break

    jobProcess = JobProcess(name = "Generating leadfield for analysis of " + str(datafile), jobProcessType = jobProcessLbexType, machine = masterMachine, processId = os.getpid(), jobStatusCode = startJobStatusCode, analysisDetail = analysisDetail,  submittedJob = submittedJob, startedTime = datetime.datetime.now()) 
    jobProcess.save()

    lbexAnalysisDetail = LBEXAnalysisDetail.objects.filter ( analysisDetail = analysisDetail ) [0]
    lbexObj = LbexObj()

    shells = Shell.objects.all()
    lbexObj.shells = shells

    lbexObj.useEGI1020Mapping = lbexAnalysisDetail.useEGI1020Mapping
    lbexObj.refLeadListChannel = lbexAnalysisDetail.refLeadListChannel

    lbexObj.scalpRadius = lbexAnalysisDetail.scalpRadius
    lbexObj.scalpConductivity = lbexAnalysisDetail.scalpConductivity
    lbexObj.sphereCenterTolerance = lbexAnalysisDetail.sphereCenterTolerance
    lbexObj.modelRadius = lbexAnalysisDetail.modelRadius
    lbexObj.modelConductivity = lbexAnalysisDetail.modelConductivity
    lbexObj.voxelSize = lbexAnalysisDetail.voxelSize
    lbexObj.coordinateOriginX = lbexAnalysisDetail.coordinateOriginX
    lbexObj.coordinateOriginY = lbexAnalysisDetail.coordinateOriginY
    lbexObj.coordinateOriginZ = lbexAnalysisDetail.coordinateOriginZ
    lbexObj.seriesCutOff = lbexAnalysisDetail.seriesCutOff
    lbexObj.tailApproxDegree = lbexAnalysisDetail.tailApproxDegree
    lbexObj.nMax = lbexAnalysisDetail.nMax
    lbexObj.concROIVolume = lbexAnalysisDetail.concROIVolume
    lbexObj.concCutOff = lbexAnalysisDetail.concCutOff
    lbexObj.concThreshold = lbexAnalysisDetail.concThreshold

    lbexFilePath = analysisDetail.lbexFilePath
    
    concEigenVectorFilePath = lbexFilePath + "/" + "cV_.txt"
    concEigenValueFilePath = lbexFilePath + "/" + "cE_.txt"

    lbexAnalysisDetail.concEigenVectorFilePath = concEigenVectorFilePath
    lbexAnalysisDetail.concEigenValueFilePath = concEigenValueFilePath

    lbexAnalysisDetail.save()

    leadList = analysisDetail.leadList
    #[ numVoxels, eegLeadField, dipoleLocations ] = createleadField(lbexObj, leadList)
    [ numVoxels, eegLeadField, dipoleLocations ] = createleadField(lbexObj, analysisDetail)

    leadFieldName = lbexAnalysisDetail.analysisDetail.leadList.headbox.code;
    lbexFilePath = lbexAnalysisDetail.analysisDetail.lbexFilePath
    
    leadFieldPath = lbexFilePath + "/" + "leadField.txt"
    dipoleLocationsPath = lbexFilePath + "/" + "dipoleLocations.txt" 
    leadFieldSVDPath = lbexFilePath + "/" + "leadFieldSVD.txt"  

    leadField = LeadField ( lbexAnalysisDetail = lbexAnalysisDetail, leadFieldName = leadFieldName, leadFieldPath = leadFieldPath, dipoleLocationsPath = dipoleLocationsPath, leadFieldSVDPath = leadFieldSVDPath, numVoxels = numVoxels ) 
    leadField.save()

    numpy.savetxt ( leadField.leadFieldPath, eegLeadField )
    numpy.savetxt ( leadField.dipoleLocationsPath, dipoleLocations )

    # leadField = LeadField.objects.get ( pk = 2 ) 

    jobProcess.jobStatusCode = endJobStatusCode
    jobProcess.completedTime = datetime.datetime.now()
    jobProcess.save()    

    # eegLeadField = numpy.genfromtxt ( leadField.leadFieldPath ) 
    # dipoleLocations = numpy.genfromtxt ( leadField.dipoleLocationsPath )

    eegLeadField = eegLeadField.transpose ()

    #print " eegLeadField[0] = " + str ( eegLeadField[0]   )
    #print " eegLeadField[1] = " + str ( eegLeadField[1]   ) 

    #print " dipoleLocations[0] = " + str ( dipoleLocations[0]   ) 
    #print " dipoleLocations[1] = " + str ( dipoleLocations[1]   ) 

    jobProcess = JobProcess(name = "Generating concentration eigenVectors and eigenValues for analysis of " + str(datafile), jobProcessType = jobProcessLbexType, machine = masterMachine, processId = os.getpid(), jobStatusCode = startJobStatusCode, analysisDetail = analysisDetail, submittedJob = submittedJob, startedTime = datetime.datetime.now()) 
    jobProcess.save()

    ddConcentration ( eegLeadField, dipoleLocations, lbexObj, lbexAnalysisDetail) 

    jobProcess.jobStatusCode = endJobStatusCode
    jobProcess.completedTime = datetime.datetime.now()
    jobProcess.save()    

    if epochType == "0": 
      updateLbexPointEvents (lbexAnalysisDetail)

    if performLbexAnalysis:
        if epochType == "0": 
            computeLbexPointEvent( analysisDetail, submittedJob, machines, startJobStatusCode, endJobStatusCode, datafile, jobProcessLbexType )
        else : 
            computeLbexAllData( analysisDetail, submittedJob, machines, startJobStatusCode, endJobStatusCode, datafile, jobProcessLbexType )

# computeSourceLocalization()

def getGridCoordinates( lbexAnalysisDetail ):

    analysisDetail = lbexAnalysisDetail.analysisDetail

    lbexLeadList = analysisDetail.leadList
    lbexHeadBox  = lbexLeadList.headbox

    useEGI1020Mapping = lbexAnalysisDetail.useEGI1020Mapping
    refLeadListChannel = lbexAnalysisDetail.refLeadListChannel

    lbexHeadBoxChannels = []

    leadListChannels = LeadListChannel.objects.filter ( leadList = lbexLeadList ) 

    lbexLeadListChannelList = []
    ten20ChannelLeadList = []

    if useEGI1020Mapping:
        ten20Channels = EGI1020Assignment.objects.all()
        for ten20Channel in ten20Channels:
            lbexLeadListChannelList.append ( ten20Channel.egiLeadListChannel ) 
            # ten20ChannelLeadList.append ( ten20Channel.tenTwentyChannel ) 

    for leadListChannel in leadListChannels:

        if useEGI1020Mapping:
            if leadListChannel not in lbexLeadListChannelList:
                continue

        EGI1020AssignmentChannel = EGI1020Assignment.objects.filter ( egiLeadListChannel = leadListChannel )[0]

        # tenTwentyChannel = EGI1020AssignmentChannel.tenTwentyChannel
        # headBoxChannel = tenTwentyChannel.headBoxChannel

        headBoxChannel = leadListChannel.headBoxChannel

        if headBoxChannel.xCoordinate == 0 and headBoxChannel.yCoordinate == 0 and headBoxChannel.zCoordinate == 0:
            continue
        sphereHeadBoxChannels.append(headBoxChannel)
        
    # refHeadboxChannels = []
    # refLeadListEGIChannel = EGI1020Assignment.objects.filter ( egiLeadListChannel = refLeadListChannel ) [0]
    # refHeadBoxChannel = refLeadListEGIChannel.tenTwentyChannel.headBoxChannel
    # refHeadboxChannels.append ( refHeadBoxChannel ) 
    refHeadboxChannels = []
    refLeadListEGIChannel = EGI1020Assignment.objects.filter ( egiLeadListChannel = refLeadListChannel ) [0]
    refHeadBoxChannel = refLeadListEGIChannel.tenTwentyChannel.headBoxChannel
    refHeadboxChannels.append ( refHeadBoxChannel ) 

    channelDistances = array ( [ sqrt(a.xCoordinate * a.xCoordinate + a.yCoordinate * a.yCoordinate + a.zCoordinate * a.zCoordinate) for a in lbexHeadBoxChannels] ) 
    sensorLocations =  [ChannelObj( a.channelId, float (lbexAnalysisDetail.scalpRadius) * ( a.xCoordinate/float(b) ), float (lbexAnalysisDetail.scalpRadius) * ( a.yCoordinate/float(b) ) , float (lbexAnalysisDetail.scalpRadius) * ( a.zCoordinate/float(b) )  )  for a,b in zip ( lbexHeadBoxChannels, channelDistances) ]

    # refChannelDistances = array ( [ sqrt(a.xCoordinate * a.xCoordinate + a.yCoordinate * a.yCoordinate + a.zCoordinate * a.zCoordinate) for a in refHeadboxChannels] ) 
    # refSensorLocations =  [ChannelObj( a.channelId, float (lbexAnalysisDetail.scalpRadius) * ( a.xCoordinate/float(b) ), float (lbexAnalysisDetail.scalpRadius) * ( a.yCoordinate/float(b) ) , float (lbexAnalysisDetail.scalpRadius) * ( a.zCoordinate/float(b) )  )  for a,b in zip ( refHeadboxChannels, channelDistances) ]
    refChannelDistances = array ( [ sqrt(a.xCoordinate * a.xCoordinate + a.yCoordinate * a.yCoordinate + a.zCoordinate * a.zCoordinate) for a in refHeadboxChannels] ) 
    refSensorLocations =  [ChannelObj( a.channelId, float (lbexAnalysisDetail.scalpRadius) * ( a.xCoordinate/float(b) ), float (lbexAnalysisDetail.scalpRadius) * ( a.yCoordinate/float(b) ) , float (lbexAnalysisDetail.scalpRadius) * ( a.zCoordinate/float(b) )  )  for a,b in zip ( refHeadboxChannels, channelDistances) ]

    # print " sensorLocations = " + str ( shape ( sensorLocations ) ) 

    xCoords = [ a.xCoordinate for a in sensorLocations ] 
    yCoords = [ a.yCoordinate for a in sensorLocations ] 
    zCoords = [ a.zCoordinate for a in sensorLocations ] 

    xLimits =  [floor( min ( xCoords ) ) , ceil ( max ( xCoords )  ) ] 
    yLimits =  [floor( min ( yCoords ) ) , ceil ( max ( yCoords )  )]  
    zLimits =  [floor( min ( zCoords ) ) , ceil ( max ( zCoords )  )]  
  
    gridLimit = [ xLimits, yLimits, zLimits ]
    fraction = GRID_LIMIT_FRACTION 

    fvs = fraction * 1

    xGrid = arange ( xLimits[0]+fvs ,   xLimits[1]-fvs, 1 )
    yGrid = arange ( yLimits[0]+fvs ,   yLimits[1]-fvs, 1 )
    zGrid = arange ( zLimits[0]+fvs ,   zLimits[1]-fvs, 1 )

    return [ xGrid, yGrid , zGrid ]

def closestSliceIndex(numberOfSections, old, user):
    indices = numpy.zeros(numberOfSections)
    indices[0] = 0
    indices[len(indices)-1] = len(old)-1
    # print " indices = " + str ( indices ) 
    for index in range (1, numberOfSections-1):
        indx = nonzero ( (old - user[index] >= 0) == 1 )[0]
        indices[index] = indx[0]
    # print " indices = " + str ( indices ) 
    
    return indices

def displaySections(xd, yd, zd, fd, numberOfSections, sectionType, parameterObj, imagePathObj, lbexFileName):

    xArray = xd.flatten(1)
    yArray = yd.flatten(1)
    zArray = zd.flatten(1)

    xmin = min(xArray) 
    xmax = max(xArray) 

    ymin = min(yArray) 
    ymax = max(yArray) 

    zmin = min(zArray) 
    zmax = max(zArray) 

    xuser = linspace(xmin, xmax, numberOfSections);
    yuser = linspace(ymin, ymax, numberOfSections);
    zuser = linspace(zmin, zmax, numberOfSections);

    # print " xuser = " + str ( xuser ) 

    Nx = len(xd ) 
    Ny = len(xd[0] ) 
    Nz = len(xd[0][0] ) 

    # print " nx = " + str ( Nx ) + " ny = " + str ( Ny ) + " nz = " + str ( Nz )  + " shape " + str ( shape ( xd))

    xold = linspace(xmin, xmax, Nx)
    yold = linspace(ymin, ymax, Ny)
    zold = linspace(zmin, zmax, Nz)

    # print " xold = " + str ( xold ) + " yold = " + str ( yold ) + " zold = " + str ( zold )  + " shape xold= " + str ( shape ( xold))
    # print " zold = " + str ( zold ) + str ( shape ( zold))

    planes = []

    if sectionType == SECTIONTYPE_TRANSVERSE:
        indices = closestSliceIndex(numberOfSections, zold, zuser)
        for index in indices:
            planes.append ( zold [ index ] )
    elif sectionType == SECTIONTYPE_CORONAL:
        indices = closestSliceIndex(numberOfSections, yold, yuser)
        for index in indices:
            planes.append ( yold [ index ] )
    elif sectionType == SECTIONTYPE_SAGITTAL:
        indices = closestSliceIndex(numberOfSections, xold, yuser)
        for index in indices:
            planes.append ( xold [ index ] )

    # print " shape fd = " + str ( shape ( fd ) ) + " fd ( 100 ) "  + str ( fd [0][:10] ) 

    fdArray = fd.flatten(1)

    # print " fdArray = " + str ( fdArray[0:20] )

    fmin = numpy.nanmin(fdArray)
    fmax = numpy.nanmax(fdArray)

    # fmin = numpy.min(fdArray)
    # fmax = numpy.max(fdArray)

    # print " fmin = " + str ( fmin ) + " fmax = " + str ( fmax ) 

    iColumns = NUMBER_COLUMNS
    iRows = ceil(numberOfSections/iColumns)
    iRows = iRows + 1
    iColumns = iColumns + 1

    # print " iColumns = " + str ( iColumns ) + " iRows = " + str ( iRows ) 

    figure()

    fontP = FontProperties()
    fontP.set_size('xx-small')

    # ax = plt.gca()
    # ax.set_axis_off()
    # lbexFileNames = os.listdir(imagePathObj.lbexPath)
    # print " spectrogramPath " + spectrogramPath

    if sectionType == SECTIONTYPE_TRANSVERSE:
        # plt.axis ( [ xmin, xmax, ymin, ymax ] )
        plt.axis('off')

        fig = plt.figure( figsize=(18, 6), dpi=100, facecolor='w', edgecolor='k')
        ax = fig.add_subplot()
        # ax.axis ( "off")

        # fig.axes.get_xaxis().set_visible(False)
        # fig.axes.get_yaxis().set_visible(False)

        lbexImageFileName = lbexFileName + "_section_transverse.png"
        lbexImageFile = open( lbexImageFileName , "w")

        #  100 volumes with 16 slices and a 32x32 inplane matrix
        # noise = N.random.randn(100,16,32,32)        

        niftiImage = [1,numberOfSections-1, xd,yd]

        for sectionIndex in range ( numberOfSections-1) :

            # print " num rows = " + str ( iRows ) + " num rows = " + str( iColumns ) + " sectionIndex = " + str( sectionIndex ) 
            plt.subplot ( int(iRows), int(iColumns), sectionIndex ) 
            sliceIndex = indices [ sectionIndex] 
            # print " section index = " + str ( sectionIndex ) + " slice index = " + str( sliceIndex ) 
            xdSection = xd[:,:,sliceIndex ]
            ydSection = yd[:,:,sliceIndex ]
            fdSection  = fd[:,:,sliceIndex ]

            niftiImage = [1,xd,yd,fd]

            plt.pcolormesh(xdSection, ydSection, fdSection, cmap='jet', vmin=fmin, vmax=fmax)
            plt.axis('off')
            # ax = plt.gca()
            # ax.set_axis_off()
        plt.savefig (lbexImageFile ,bbox_inches='tight') 

    elif sectionType == SECTIONTYPE_CORONAL:
        # plt.axis ( [ xmin, xmax, zmin, zmax ] )
        plt.axis('off')
        
        fig = plt.figure( figsize=(18, 6), dpi=100, facecolor='w', edgecolor='k')
        ax = fig.add_subplot()
        # ax.axis ( "off")
        # fig.axes.get_xaxis().set_visible(False)
        # fig.axes.get_yaxis().set_visible(False)

        lbexImageFileName = lbexFileName + "_section_coronal.png"
        lbexImageFile = open( lbexImageFileName , "w")
        
        for sectionIndex in range ( numberOfSections-1) :
            # print " num rows = " + str ( iRows ) + " num rows = " + str( iColumns ) + " sectionIndex = " + str( sectionIndex ) 
            plt.subplot ( int(iRows), int(iColumns), sectionIndex ) 
            sliceIndex = indices [ sectionIndex] 
            # print " section index = " + str ( sectionIndex ) + " slice index = " + str( sliceIndex ) 
            xdSection = xd[sliceIndex,:,: ]
            zdSection = zd[sliceIndex,:,: ]
            fdSection  = fd[sliceIndex,:,: ]

            niftiImage = [1,xd,yd,fd]

            plt.pcolormesh(xdSection, zdSection, fdSection, cmap='jet', vmin=fmin, vmax=fmax)
            plt.axis('off')
            # ax = plt.gca()
            # ax.set_axis_off()
        plt.savefig ( lbexImageFile ,bbox_inches='tight' ) 

    elif sectionType == SECTIONTYPE_SAGITTAL:

        plt.axis ( [ ymin, ymax, zmin, zmax ] )
        plt.axis('off')
        
        fig = plt.figure( figsize=(18, 6), dpi=100, facecolor='w', edgecolor='k')
        ax = fig.add_subplot()
        # ax.axis ( "off")

        # fig.axes.get_xaxis().set_visible(False)
        # fig.axes.get_yaxis().set_visible(False)

        lbexImageFileName = lbexFileName + "_section_sagittal.png"
        lbexImageFile = open( lbexImageFileName , "w")

        for sectionIndex in range ( numberOfSections-1) :
            # print " num rows = " + str ( iRows ) + " num rows = " + str( iColumns ) + " sectionIndex = " + str( sectionIndex ) 
            plt.subplot ( int(iRows), int(iColumns), sectionIndex ) 
            sliceIndex = indices [ sectionIndex] 
            # print " section index = " + str ( sectionIndex ) + " slice index = " + str( sliceIndex ) 
            ydSection = yd[:,sliceIndex,: ]
            zdSection = zd[:,sliceIndex, : ]
            fdSection  = fd[:,sliceIndex,: ]

            niftiImage = [1,xd,yd,fd]

            plt.pcolormesh(ydSection, zdSection, fdSection, cmap='jet', vmin=fmin, vmax=fmax)
            plt.axis('off')
        plt.savefig ( lbexImageFile ,bbox_inches='tight' ) 

# class VolumeSlicer(HasTraits):
#     # The data to plot
#     data = Array()

#     # The 4 views displayed
#     scene3d = Instance(MlabSceneModel, ())
#     scene_x = Instance(MlabSceneModel, ())
#     scene_y = Instance(MlabSceneModel, ())
#     scene_z = Instance(MlabSceneModel, ())

#     # The data source
#     data_src3d = Instance(Source)

#     # The image plane widgets of the 3D scene
#     ipw_3d_x = Instance(PipelineBase)
#     ipw_3d_y = Instance(PipelineBase)
#     ipw_3d_z = Instance(PipelineBase)

#     _axis_names = dict(x=0, y=1, z=2)


#     #---------------------------------------------------------------------------
#     def __init__(self, **traits):
#         super(VolumeSlicer, self).__init__(**traits)
#         # Force the creation of the image_plane_widgets:
#         self.ipw_3d_x
#         self.ipw_3d_y
#         self.ipw_3d_z


#     #---------------------------------------------------------------------------
#     # Default values
#     #---------------------------------------------------------------------------
#     def _data_src3d_default(self):
#         return mlab.pipeline.scalar_field(self.data,
#                             figure=self.scene3d.mayavi_scene)

#     def make_ipw_3d(self, axis_name):
#         ipw = mlab.pipeline.image_plane_widget(self.data_src3d,
#                         figure=self.scene3d.mayavi_scene,
#                         plane_orientation='%s_axes' % axis_name)
#         return ipw

#     def _ipw_3d_x_default(self):
#         return self.make_ipw_3d('x')

#     def _ipw_3d_y_default(self):
#         return self.make_ipw_3d('y')

#     def _ipw_3d_z_default(self):
#         return self.make_ipw_3d('z')


#     #---------------------------------------------------------------------------
#     # Scene activation callbaks
#     #---------------------------------------------------------------------------
#     @on_trait_change('scene3d.activated')
#     def display_scene3d(self):
#         outline = mlab.pipeline.outline(self.data_src3d,
#                         figure=self.scene3d.mayavi_scene,
#                         )
#         self.scene3d.mlab.view(40, 50)
#         # Interaction properties can only be changed after the scene
#         # has been created, and thus the interactor exists
#         for ipw in (self.ipw_3d_x, self.ipw_3d_y, self.ipw_3d_z):
#             # Turn the interaction off
#             ipw.ipw.interaction = 0
#         self.scene3d.scene.background = (0, 0, 0)
#         # Keep the view always pointing up
#         self.scene3d.scene.interactor.interactor_style = \
#                                  tvtk.InteractorStyleTerrain()


#     def make_side_view(self, axis_name):
#         scene = getattr(self, 'scene_%s' % axis_name)

#         # To avoid copying the data, we take a reference to the
#         # raw VTK dataset, and pass it on to mlab. Mlab will create
#         # a Mayavi source from the VTK without copying it.
#         # We have to specify the figure so that the data gets
#         # added on the figure we are interested in.
#         outline = mlab.pipeline.outline(
#                             self.data_src3d.mlab_source.dataset,
#                             figure=scene.mayavi_scene,
#                             )
#         ipw = mlab.pipeline.image_plane_widget(
#                             outline,
#                             plane_orientation='%s_axes' % axis_name)
#         setattr(self, 'ipw_%s' % axis_name, ipw)

#         # Synchronize positions between the corresponding image plane
#         # widgets on different views.
#         ipw.ipw.sync_trait('slice_position',
#                             getattr(self, 'ipw_3d_%s'% axis_name).ipw)

#         # Make left-clicking create a crosshair
#         ipw.ipw.left_button_action = 0
#         # Add a callback on the image plane widget interaction to
#         # move the others
#         def move_view(obj, evt):
#             position = obj.GetCurrentCursorPosition()
#             for other_axis, axis_number in self._axis_names.iteritems():
#                 if other_axis == axis_name:
#                     continue
#                 ipw3d = getattr(self, 'ipw_3d_%s' % other_axis)
#                 ipw3d.ipw.slice_position = position[axis_number]

#         ipw.ipw.add_observer('InteractionEvent', move_view)
#         ipw.ipw.add_observer('StartInteractionEvent', move_view)

#         # Center the image plane widget
#         ipw.ipw.slice_position = 0.5*self.data.shape[
#                     self._axis_names[axis_name]]

#         # Position the view for the scene
#         views = dict(x=( 0, 90),
#                      y=(90, 90),
#                      z=( 0,  0),
#                      )
#         scene.mlab.view(*views[axis_name])
#         # 2D interaction: only pan and zoom
#         scene.scene.interactor.interactor_style = \
#                                  tvtk.InteractorStyleImage()
#         scene.scene.background = (0, 0, 0)


#     @on_trait_change('scene_x.activated')
#     def display_scene_x(self):
#         return self.make_side_view('x')

#     @on_trait_change('scene_y.activated')
#     def display_scene_y(self):
#         return self.make_side_view('y')

#     @on_trait_change('scene_z.activated')
#     def display_scene_z(self):
#         return self.make_side_view('z')


#     #---------------------------------------------------------------------------
#     # The layout of the dialog created
#     #---------------------------------------------------------------------------
#     view = View(HGroup(
#                   Group(
#                        Item('scene_y',
#                             editor=SceneEditor(scene_class=Scene),
#                             height=250, width=300),
#                        Item('scene_z',
#                             editor=SceneEditor(scene_class=Scene),
#                             height=250, width=300),
#                        show_labels=False,
#                   ),
#                   Group(
#                        Item('scene_x',
#                             editor=SceneEditor(scene_class=Scene),
#                             height=250, width=300),
#                        Item('scene3d',
#                             editor=SceneEditor(scene_class=MayaviScene),
#                             height=250, width=300),
#                        show_labels=False,
#                   ),
#                 ),
#                 resizable=True,
#                 title='Volume Slicer',
#                 )

def displayData(data, lbexAnalysisDetail, imagePathObj, lbexFileName, lbexSelectionObj):

 try:
     
     [xGrid  , yGrid , zGrid] = getGridCoordinates(lbexAnalysisDetail)

     # print " len x = " + str ( len(xArray) ) + " len y = " + str ( len(yArray) )  + " len   z = " +  str ( len(zArray) ) 
    
     # print " x = " + str ( xArray[0:25] ) + "  y = " + str ( yArray[0:25]  )  + " z = " +  str ( zArray[0:25]  ) 

     nCubes = len(xGrid) * len(yGrid) * len(zGrid)

     x = zeros(nCubes)
     y = zeros(nCubes)
     z = zeros(nCubes)
     
     inside = zeros(nCubes)

     xArray = []
     yArray = []
     zArray = []

     [x,y,z] = meshgrid( xGrid  , yGrid , zGrid )
     # print " shape   = " + str ( shape ( x  ) ) + " - " +  str ( shape ( y  ) ) + " - " + str ( shape ( z  ) )
     xArray = x.flatten(1)
     yArray = y.flatten(1)
     zArray = z.flatten(1)
     # print " max max   = " + str ( max ( xArray  ) ) + " - " +  str ( max ( yArray  ) ) + " - " + str ( max ( zArray  ) )
     radialDistanceArray = [ a*a + b*b + c*c for a,b,c in zip ( xArray, yArray, zArray ) ]

     inside = [index for index,elem in enumerate(radialDistanceArray) if elem <= lbexAnalysisDetail.scalpRadius * lbexAnalysisDetail.scalpRadius]

     # print " inside = " + str ( len ( inside ) ) 
     # print " nCubes = " + str ( nCubes ) 

     xd = reshape(x, [len(yGrid),len(xGrid),len(zGrid)])
     yd = reshape(y, [len(yGrid),len(xGrid),len(zGrid)])
     zd = reshape(z, [len(yGrid),len(xGrid),len(zGrid)])

     # print " xd shape  = " + str ( shape ( xd  ) ) 

     fd = numpy.zeros(nCubes)
     fd.fill(numpy.nan)     

     # fd = [ y for x,y in enumerate ( zip ( data, inside ) ) if x == y[1]  ]
     for index, indexVal in enumerate ( inside )  :
         # print " index = " + str ( index ) + " type = " + str ( type ( data  )  ) + " value = " + str ( data [index]  ) 
         fd[indexVal] = data[index]

     fd =  reshape ( fd, [len(yGrid), len(xGrid), len(zGrid)],order="F") 
     # print " fd  = " + str ( fd  [:2,:11,0] ) 
     
     parameterObj = ParameterObj ( OPACITY_INCREASE, OPACITY_DECREASE, COLOR_MAP ) 

     show2DSections = lbexSelectionObj.show2DSections
     show2DTransverseView = lbexSelectionObj.show2DTransverseView
     show2DCoronalView = lbexSelectionObj.show2DCoronalView
     show2DSagittalView = lbexSelectionObj.show2DSagittalView

     numLbex2DSections = lbexSelectionObj.numLbex2DSections

     show3DSections = lbexSelectionObj.show3DSections

     if show2DSections:

         if show2DTransverseView:

             sectionType = SECTIONTYPE_TRANSVERSE 
             displaySections( xd, yd, zd, fd, numLbex2DSections, sectionType, parameterObj, imagePathObj, lbexFileName )

         if show2DCoronalView:

             sectionType = SECTIONTYPE_CORONAL 
             displaySections( xd, yd, zd, fd, numLbex2DSections, sectionType, parameterObj, imagePathObj, lbexFileName )

         if show2DSagittalView:

             sectionType = SECTIONTYPE_SAGITTAL
             displaySections( xd, yd, zd, fd, numLbex2DSections, sectionType, parameterObj, imagePathObj, lbexFileName )

     fig = plt.figure()

     if show3DSections:

         m = VolumeSlicer(data=data)
         m.configure_traits()

 except:
     traceback.print_exc(file=sys.stdout)
 return

def generateLbexImages(lbexAnalysisDetail, imagePathObj, plt, lbexSelectionObjs ):

    try:

        analysisDetail = lbexAnalysisDetail.analysisDetail
        datafile = analysisDetail.datafile
        lbexImagePath = analysisDetail.lbexImagePath

        for lbexSelectionObj in lbexSelectionObjs:
            
            lbexDataFileName = lbexImagePath + "/" + datafile.name + "_" + str ( lbexSelectionObj.lbexPointEventId ) 

            print (" #### lbexDataFileName = " + lbexDataFileName)
            data = []

            f = open ( lbexDataFileName ,"r")
           
            for line in f:
                ar = line.replace("\n", "").split(",") 

                data.extend( ar )
                # print " *************** for index " + str ( index ) + " data = " + str ( type (data)  ) + " shape = " + str ( shape ( data ) )
            
            print (" in display data LBEX  " )

            displayData( array(data) , lbexAnalysisDetail, imagePathObj, lbexDataFileName, lbexSelectionObj)
    except:
     
        traceback.print_exc(file=sys.stdout)

    return
       
# displayData() 

# data = []
# f = open ("/home/hduser/cornell/data/lbex_data/data.txt_76","r")
# for line in f:
#     ar = line.split(",")
#     for a in ar:
#         if a!= "":
#             data.append(float(a))
# print " len data " + str ( len ( data ) )
# data = array(data)

# analysisDetail = AnalysisDetail.objects.get ( pk = 45)
# headbox = analysisDetail.leadList.headbox     
# lbexAnalysisDetail = LBEXAnalysisDetail.objects.get ( pk = 26)
# xGrid  , yGrid , zGrid = getGridCoordinates(headbox, lbexAnalysisDetail)
# print " type data = " + str ( type ( data )) + " type x = " + str ( type ( xGrid )) + " size data = " + str (  data.shape  )  

# nCubes = len(xGrid) * len(yGrid) * len(zGrid)

# x = zeros(nCubes)
# y = zeros(nCubes)
# z = zeros(nCubes)
     
# inside = zeros(nCubes)

# xArray = []
# yArray = []
# zArray = []

# [xg,yg,zg] = meshgrid( xGrid  , yGrid , zGrid )

# xArray = xg.flatten(1)
# yArray = yg.flatten(1)
# zArray = zg.flatten(1)

# radialDistanceArray = [ a*a + b*b + c*c for a,b,c in zip ( xArray, yArray, zArray ) ]

# inside = [index for index,elem in enumerate(radialDistanceArray) if elem<=SCALP_RADIUS * SCALP_RADIUS]

# xd = reshape(x, [len(yGrid),len(xGrid),len(zGrid)])
# yd = reshape(y, [len(yGrid),len(xGrid),len(zGrid)])
# zd = reshape(z, [len(yGrid),len(xGrid),len(zGrid)])

# fd = np.zeros(nCubes)
# # fd.fill(numpy.nan)     

# # fd = [ y for x,y in enumerate ( zip ( data, inside ) ) if x == y[1]  ]
# for index, indexVal in enumerate ( inside )  :
#     print " index = " + str ( index ) + " index val " + str ( indexVal ) 
#     if index > ( len ( data ) -1 ) :
#        break
#     fd[indexVal] = data[index]

# data =  reshape ( fd, [len(yGrid), len(xGrid), len(zGrid)],order="F") 
     # print " fd  = " + str ( fd  [:2,:11,0] ) 
################################################################################
# The object implementing the dialog


# def displaySections3D(xd, yd, zd, fd, numberOfSections, sectionType, parameterObj):

#     xArray1stDim = len(xd)
#     xArray2ndDim = len(xd[0])
#     xArray3rdDim = len(xd[0][0])

#     xArray = xd.flatten(1)
#     yArray = yd.flatten(1)
#     zArray = zd.flatten(1)
#     fArray = fd.flatten(1)

#     xmin = min(xArray) 
#     xmax = max(xArray) 

#     ymin = min(yArray) 
#     ymax = max(yArray) 

#     zmin = min(zArray) 
#     zmax = max(zArray) 

#     xuser = linspace(xmin, xmax, numberOfSections);
#     yuser = linspace(ymin, ymax, numberOfSections);
#     zuser = linspace(zmin, zmax, numberOfSections);

#     fontP = FontProperties()
#     fontP.set_size('xx-small')

#     ax = plt.gca(projection='3d')
#     # ax.set_axis_off()

#     # ax = p3.Axes3D(fig)

#     if sectionType == SECTIONTYPE_TRANSVERSE:
#         numberOfSections = min( numberOfSections, xArray3rdDim )
#         zplanes = linspace( zmin, zmax, numberOfSections)
#         # print " zplanes = " + str ( zplanes ) 
#         print " shape x = " + str ( shape ( xArray ) ) + " shape y = " + str ( shape ( yArray ) )  + " shape f = " + str ( shape ( fArray ) ) 
#         print " x = " + str ( xArray[:10] ) + " y = " + str ( yArray[:10] ) + " shape f = " + str (  fArray[:50] )
#         ax.plot_surface(xArray, yArray, fArray)
#         ax.set_zlim(zmin, zmax)
#         plt.show()
#         # ax.contourf3D(X,Y,Z)
#         # hs = slice(xd,yd,zd,fd, [], [], zplanes );
#         # title = 'Transverse'
#     elif sectionType == SECTIONTYPE_CORONAL:
#         numberOfSections = min( numberOfSections, xArray1stDim )
#         yplanes = linspace( ymin, ymax, numberOfSections)
#         # hs = slice(xd,yd,zd,fd, [], [], zplanes );
#         # title = 'Coronal'
#     elif sectionType == SECTIONTYPE_TRANSVERSE:
#        numberOfSections = min( numberOfSections, xArray2ndDim )
#        xplanes = linspace( xmin, xmax, numberOfSections)
#         # hs = slice(xd,yd,zd,fd, [], [], zplanes );
#        # title = 'Sagittal'
