# This script is provided by Adam Deller

import sys, math, os

####### FIND A LINE STARTING WITH A GIVEN KEY ##################################
def getToLine(key, imlines, linenum):
    while linenum < len(imlines) and not imlines[linenum].split(':')[0] == key:
        linenum += 1
    if linenum == len(imlines):
        print "Couldn't find key " + key + " - aborting!"
        sys.exit(1)
    return linenum, imlines[linenum].split(':')[1]

####### CONVERT A YEAR MONTH DAY DATE TO MJD ###################################
def ymd2mjd(year, month, day):
    return year*367 - int(7*(year + int((month + 9)/12))/4) + \
           int(275*month/9) + day - 678987


################################################################################
# Container classes
################################################################################
class VexScan:
    def __init__(self, scanname, starttime, stoptime, source):
        self.scanname  = scanname
        self.starttime = starttime
        self.stoptime  = stoptime
        self.source    = source
        self.startmjd  = getVexMJD(starttime)
        self.stopmjd   = getVexMJD(stoptime)

    def incStopTime(self, incsecs):
        self.stopmjd += incsecs/86400.0

    def setStopTime(self, stoptime):
        self.stoptime = stoptime
        self.stopmjd = getVexMJD(stoptime)

class Telescope:
    def __init__(self, name, x, y, z, axisoff):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.axisoff = axisoff

class Source:
    def __init__(self, name, ra, dec):
        self.name = name
        self.ra = ra
        self.dec = dec

    def getSepArcmin(self, rarad, decrad):
        sinsqdecdiff = math.sin((decrad-self.dec)/2.0)
        sinsqdecdiff = sinsqdecdiff*sinsqdecdiff
        sinsqradiff  = math.sin((rarad-self.ra)/2.0)
        sinsqradiff  = sinsqradiff*sinsqradiff

        return 180.0*60.0*2*math.asin(math.sqrt(sinsqdecdiff +
                           math.cos(decrad)*math.cos(self.dec)*sinsqradiff))/math.pi

class Scan:
    def __init__(self, imlines, calclines, ilinenum, clinenum, numantennas,
                 polyorder, polyinterval, isdifx15, sources, startmjd, startsec,
                 compatibilityMode):
        self.startmjd = startmjd
        self.startsec = startsec
        self.polyorder = polyorder
        self.polyinterval = polyinterval
        self.polymjd = []
        self.polysec = []
        self.delays  = []
        self.dry  = []
        self.wet  = []
        self.u       = []
        self.v       = []
        self.w       = []
        self.source = []

        if compatibilityMode: #"im" is really delay, "calc" is really uvw.
                              #Initialise accordingly
            self.initialiseCompatibility(imlines, calclines, ilinenum, clinenum,
                                         numantennas, polyorder, polyinterval,
                                         sources, startmjd, startsec)
            return

        # Do the calc side first
        self.scanstart = int(calclines[clinenum+1].split(':')[1].strip())
        if isdifx15:
            self.scandur = int(calclines[clinenum+0].split(':')[1].strip())
            name = calclines[clinenum+3].split(':')[1].strip()
            ra = float(calclines[clinenum+4].split(':')[1])
            dec = float(calclines[clinenum+5].split(':')[1])
            self.source.append(Source(name, ra, dec))
            clinenum += 8
            self.numphasecentres = 0
        else:
            self.scandur = int(calclines[clinenum+2].split(':')[1])
            self.numphasecentres = int(calclines[clinenum+7].split(':')[1])
            self.source.append(sources[int(calclines[clinenum+6].split(':')[1])])
            for i in range(self.numphasecentres):
                self.source.append(sources[int(calclines[clinenum+8+i].split(':')[1])])
            clinenum += self.numphasecentres+8

        self.startsec += self.scanstart

        # Then then im side
        if isdifx15:
            self.numpoly = int(imlines[ilinenum+1].split(':')[1])
            ilinenum += 2
        else:
            self.numpoly = int(imlines[ilinenum+self.numphasecentres+2].split(':')[1])
            ilinenum += self.numphasecentres + 3
        for i in range(self.numpoly):
            self.polymjd.append(int(imlines[ilinenum].split(':')[1]))
            ilinenum += 1
            self.polysec.append(int(imlines[ilinenum].split(':')[1]))
            ilinenum += 1
            self.delays.append([])
            self.dry.append([])
            self.wet.append([])
            self.u.append([])
            self.v.append([])
            self.w.append([])
            for j in range(self.numphasecentres+1):
                self.delays[i].append([])
                self.dry[i].append([])
                self.wet[i].append([])
                self.u[i].append([])
                self.v[i].append([])
                self.w[i].append([])
                for k in range(numantennas):
                    self.delays[i][j].append([])
                    self.dry[i][j].append([])
                    self.wet[i][j].append([])
                    self.u[i][j].append([])
                    self.v[i][j].append([])
                    self.w[i][j].append([])
                    d = (imlines[ilinenum].split(':')[1]).split()
                    r = (imlines[ilinenum+1].split(':')[1]).split()
                    t = (imlines[ilinenum+2].split(':')[1]).split()
                    u = (imlines[ilinenum+5].split(':')[1]).split()
                    v = (imlines[ilinenum+6].split(':')[1]).split()
                    w = (imlines[ilinenum+7].split(':')[1]).split()
                    for l in range(polyorder+1):
                        self.delays[i][j][k].append(float(d[l]))
                        self.dry[i][j][k].append(float(r[l]))
                        self.wet[i][j][k].append(float(t[l]))
                        self.u[i][j][k].append(float(u[l]))
                        self.v[i][j][k].append(float(v[l]))
                        self.w[i][j][k].append(float(w[l]))
                    ilinenum += 8

        # Finally store the end line numbers
        self.endclinenum = clinenum
        self.endilinenum = ilinenum

    def initialiseCompatibility(self, delaylines, uvwlines, dlinenum, ulinenum,
                                numantennas, polyorder, polyinterval,
                                sources, startmjd, startsec):
        self.scanstart = int(delaylines[dlinenum+1].split(':')[1].strip())
        self.scandur   = int(delaylines[dlinenum+0].split(':')[1].strip())
        name = uvwlines[ulinenum+2].split(':')[1].strip()
        ra   = float(uvwlines[ulinenum+3].split(':')[1])
        dec  = float(uvwlines[ulinenum+4].split(':')[1])
        self.source.append(Source(name, ra, dec))
        dlinenum += 4
        ulinenum += 6
        self.numphasecentres = 0
        self.startsec += self.scanstart
        self.numpoly = self.scandur

        for i in range(self.numpoly):
            lastdelays    = delaylines[i + dlinenum - 1][20:].split()
            currentdelays = delaylines[i + dlinenum][20:].split()
            nextdelays    = delaylines[i + dlinenum + 1][20:].split()
            lastuvws      = uvwlines[i + ulinenum - 1][20:].split()
            currentuvws   = uvwlines[i + ulinenum][20:].split()
            nextuvws      = uvwlines[i + ulinenum + 1][20:].split()
            self.polymjd.append(self.startmjd)
            self.polysec.append(self.startsec + i)
            self.delays.append([])
            self.u.append([])
            self.v.append([])
            self.w.append([])
            for j in range(self.numphasecentres+1):
                self.delays[i].append([])
                self.u[i].append([])
                self.v[i].append([])
                self.w[i].append([])
                for k in range(numantennas):
                    d0 = float(lastdelays[k])
                    d1 = float(currentdelays[k])
                    d2 = float(nextdelays[k])
                    u0 = float(lastuvws[k*3+0])
                    u1 = float(currentuvws[k*3+0])
                    u2 = float(nextuvws[k*3+0])
                    v0 = float(lastuvws[k*3+1])
                    v1 = float(currentuvws[k*3+1])
                    v2 = float(nextuvws[k*3+1])
                    w0 = float(lastuvws[k*3+2])
                    w1 = float(currentuvws[k*3+2])
                    w2 = float(nextuvws[k*3+2])
                    self.delays[i][j].append([])
                    self.u[i][j].append([])
                    self.v[i][j].append([])
                    self.w[i][j].append([])
                    self.delays[i][j][k].append(d1)
                    self.delays[i][j][k].append(d2/2.0 - d0/2.0)
                    self.delays[i][j][k].append(d2/2.0 + d0/2.0 - d1)
                    #if i == 0 or i == self.numpoly-1:
                    #    print self.delays[i][j][k]
                    self.u[i][j][k].append(u1)
                    self.u[i][j][k].append(u2/2.0 - u0/2.0)
                    self.u[i][j][k].append(u2/2.0 + u0/2.0 - u1)
                    self.v[i][j][k].append(v1)
                    self.v[i][j][k].append(v2/2.0 - v0/2.0)
                    self.v[i][j][k].append(v2/2.0 + v0/2.0 - v1)
                    self.w[i][j][k].append(w1)
                    self.w[i][j][k].append(w2/2.0 - w0/2.0)
                    self.w[i][j][k].append(w2/2.0 + w0/2.0 - w1)

        # Finally store the end line numbers
        self.enddlinenum = dlinenum + self.scandur + 2
        self.endulinenum = ulinenum + self.scandur + 2

    def getSourceIndex(self, sourcename):
        for i in range(self.numphasecentres+1):
            #print "Comparing *" + self.source[i].name.strip() + "* with *" + sourcename.strip() + "*"
            #print self.source[i].name.strip() == sourcename.strip()
            if self.source[i].name.strip() == sourcename.strip():
                return i
        if self.source[i].name.strip() == "SCAN_GAP":
            return -1
        else:
            print "Couldn't find source " + sourcename.strip() + " - my pointing centre was " + self.source[0].name.strip()
            return -999

    def containsTime(self, mjd, sec):
#        print "my mjd and sec is %d, %d, being asked for %d, %f" % (self.startmjd, self.startsec, mjd, sec)
        offset = (mjd-self.startmjd)*86400 + sec - self.startsec
        if offset >= 0 and offset < self.scandur + 0.1: #Crappy rpfits time stamps!
            return True
        return False

    def getPolyEntryAndOffset(self, mjd, second):
        polyentry = -1
        for i in range(self.numpoly):
            polyoffset = (mjd - self.polymjd[i])*86400 + second - self.polysec[i]
            #print "mjd: " + str(mjd), "polymjd: " + str(self.polymjd[i]), "second: " + str(second), "polysec: " + str(self.polysec[i]), "polyoffset: " + str(polyoffset) 
            if polyoffset >= 0 and polyoffset <= self.polyinterval:
                polyentry = i
                print polyentry
                break
        if polyentry < 0:
            # print polyoffset
            print "Couldn't find a poly entry for time " + str(mjd) + ", sec " + \
                  str(second) + " in this scan - aborting!"
            polyoffset = 0
            
        return polyentry, polyoffset

    def getUVW(self, mjd, second, antenna1index, antenna2index, sourceindex):
        toreturn = [0.0, 0.0, 0.0]
        if sourceindex > self.numphasecentres:
            print "Trying to get source index " + str(sourceindex) + " when this " + \
                  "scan only has " + str(self.numphasecentres) + " phase centres"
            return toreturn
        polyentry, polyoffset = self.getPolyEntryAndOffset(mjd, second)
        if polyentry < 0:
            print "Couldn't find a poly entry for time " + str(mjd) + ", sec " + \
                  str(second) + " in this scan"
            return toreturn
        xval = 1.0
        for i in range(self.polyorder+1):
            toreturn[0] += self.u[polyentry][sourceindex][antenna1index][i]*xval
            toreturn[1] += self.v[polyentry][sourceindex][antenna1index][i]*xval
            toreturn[2] += self.w[polyentry][sourceindex][antenna1index][i]*xval
            xval *= polyoffset
        xval = 1.0
        for i in range(self.polyorder+1):
            toreturn[0] -= self.u[polyentry][sourceindex][antenna2index][i]*xval
            toreturn[1] -= self.v[polyentry][sourceindex][antenna2index][i]*xval
            toreturn[2] -= self.w[polyentry][sourceindex][antenna2index][i]*xval
            xval *= polyoffset
        return toreturn

    def getDelayAndRate(self, mjd, second, antennaindex, sourceindex):
        if sourceindex > self.numphasecentres:
            print "Trying to get source index " + str(sourceindex) + " when this " + \
                  "scan only has " + str(self.numphasecentres) + " phase centres"
            return 0.0, 0.0
        polyentry, polyoffset = self.getPolyEntryAndOffset(mjd, second)
        if polyentry < 0:
            print "Couldn't find a poly entry for time " + str(mjd) + ", sec " + \
                  str(second) + " in this scan"
            return 0.0, 0.0
        correction = 0.0
        delay = 0.0
        rate = 0.0
        accel = 0.0
        xval = 1.0
        #print "Polyentry is " + str(polyentry) + ", polyoffset is " + str(polyoffset)
        #print self.delays[polyentry][sourceindex][antennaindex]
        #print self.dry[polyentry][sourceindex][antennaindex]
        #print self.wet[polyentry][sourceindex][antennaindex]
        for i in range(self.polyorder+1):
            delay += self.delays[polyentry][sourceindex][antennaindex][i]*xval
            #delay += self.dry[polyentry][sourceindex][antennaindex][i]*xval
            #delay += self.wet[polyentry][sourceindex][antennaindex][i]*xval
            if polyoffset > 0.0:
              rate += i*self.delays[polyentry][sourceindex][antennaindex][i]*xval/polyoffset
              accel+= 0.5*i*(i-1)*self.delays[polyentry][sourceindex][antennaindex][i]*xval/(polyoffset**2)
              #rate += i*self.dry[polyentry][sourceindex][antennaindex][i]*xval/polyoffset
              #rate += i*self.wet[polyentry][sourceindex][antennaindex][i]*xval/polyoffset
            else:
              rate = self.delays[polyentry][sourceindex][antennaindex][1]
              accel = self.delays[polyentry][sourceindex][antennaindex][2]
              #rate += self.dry[polyentry][sourceindex][antennaindex][1]
              #rate += self.wet[polyentry][sourceindex][antennaindex][1]
            xval *= polyoffset
        #delay += delay*rate/1e6
        #print "Calculated delay is %10.6f, calculated rate is %10.6f, xval is %10.6f" % (delay, rate, xval)
        #return delay/1000000.0, rate/1000000.0 #microsec -> seconds and us/sec -> sec/sec
        return delay, rate, accel #microsec, us/sec, and us/sec^2

class Model:
    def __init__(self, calcfile, delayfile=None, uvwfile=None):
        if calcfile == None:
            if not os.path.exists(delayfile):
                print "Delay file (compatibility mode) does not exist- aborting!"
                sys.exit()
            if not os.path.exists(uvwfile):
                print "UVW file (compatibility mode) does not exist- aborting!"
                sys.exit()
            self.initialiseCompatibility(delayfile, uvwfile)
            return
        if not os.path.exists(calcfile):
            print "Calc file " + calcfile + " does not exist - aborting"
            sys.exit()
        self.compatibilityMode = False
        calcin = open(calcfile)
        self.filename = calcfile
        calclines = calcin.readlines()
        calcin.close()
        clinenum = 0
        clinenum, val = getToLine("DIFX VERSION", calclines, clinenum)
        difxverstr = val.strip()
        clinenum, val = getToLine("IM FILENAME", calclines, clinenum)
        imfilename = val.strip()
        self.isDifx15 = False
        if "1.5" in difxverstr:
            self.isDifx15 = True
        imin = open(imfilename)
        imlines = imin.readlines()
        imin.close()
        ilinenum = 0
        ilinenum, yearstr = getToLine("START YEAR", imlines, ilinenum)
        ilinenum, monthstr = getToLine("START MONTH", imlines, ilinenum)
        ilinenum, daystr = getToLine("START DAY", imlines, ilinenum)
        ilinenum, hourstr = getToLine("START HOUR", imlines, ilinenum)
        ilinenum, minstr = getToLine("START MINUTE", imlines, ilinenum)
        ilinenum, secstr = getToLine("START SECOND", imlines, ilinenum)
        self.startmjd = ymd2mjd(int(yearstr), int(monthstr), int(daystr))
        self.startsec = int(hourstr)*3600 + int(minstr)*60 + int(secstr)
        ilinenum, polyorderstr = getToLine("POLYNOMIAL ORDER", imlines, ilinenum)
        self.polyorder = int(polyorderstr)
        ilinenum, polyintervalstr = getToLine("INTERVAL (SECS)", imlines, ilinenum)
        self.polyinterval = int(polyintervalstr)
        ilinenum, numtelstr = getToLine("NUM TELESCOPES", imlines, ilinenum)
        self.numantennas = int(numtelstr)
        clinenum = 0
        clinenum, numtelstr = getToLine("NUM TELESCOPES", calclines, clinenum)
        if self.numantennas != int(numtelstr):
            print "Number of telescopes doesn't match between calc and IM files!"
            sys.exit()
        self.antennas = []
        self.antennamap = {}
        for i in range(self.numantennas):
            clinenum, val = getToLine("TELESCOPE %d NAME" % i, calclines, clinenum)
            antennaname = val.strip()
            clinenum, val = getToLine("TELESCOPE %d OFFSET (m)" % i, calclines, clinenum)
            axisoffset = float(val)
            clinenum, val = getToLine("TELESCOPE %d X (m)" % i, calclines, clinenum)
            x = float(val.strip())
            clinenum, val = getToLine("TELESCOPE %d Y (m)" % i, calclines, clinenum)
            y = float(val.strip())
            clinenum, val = getToLine("TELESCOPE %d Z (m)" % i, calclines, clinenum)
            z = float(val.strip())
            self.antennas.append(Telescope(antennaname, x, y, z, axisoffset))
            self.antennamap[antennaname] = i
        self.sources = []
        if not self.isDifx15:
            clinenum, val = getToLine("NUM SOURCES", calclines, clinenum)
            self.numsources = int(val)
            for i in range(self.numsources):
                clinenum, val = getToLine("SOURCE %d NAME" % i, calclines, clinenum)
                name = val.strip()
                clinenum, val  = getToLine("SOURCE %d RA" % i, calclines, clinenum)
                ra = float(val.strip())
                clinenum, val  = getToLine("SOURCE %d DEC" % i, calclines, clinenum)
                dec = float(val.strip())
                self.sources.append(Source(name, ra, dec))
        clinenum, val = getToLine("NUM SCANS", calclines, clinenum)
        self.numscans = int(val)
        ilinenum, val = getToLine("NUM SCANS", imlines, ilinenum)
        if not self.numscans == int(val):
            print "Error - .im and .calc file disagree on number of scans!"
            sys.exit(1)
        clinenum += 1
        ilinenum += 1
        self.scans = []
        for j in range(self.numscans):
            toadd = Scan(imlines, calclines, ilinenum, clinenum,
                         self.numantennas, self.polyorder, self.polyinterval,
                         self.isDifx15, self.sources, self.startmjd, self.startsec,
                         self.compatibilityMode)
            self.scans.append(toadd)
            ilinenum = toadd.endilinenum
            clinenum = toadd.endclinenum
        self.lastscannum = 0
        self.changed = True

    def initialiseCompatibility(self, delayfile, uvwfile):
        self.isDifx15 = False
        self.compatibilityMode = True
        delayin = open(delayfile)
        delaylines = delayin.readlines()
        delayin.close()
        uvwin = open(uvwfile)
        uvwlines = uvwin.readlines()
        uvwin.close()
        dlinenum = 0
        ulinenum = 0
        dlinenum, yearstr = getToLine("START YEAR", delaylines, dlinenum)
        dlinenum, monthstr = getToLine("START MONTH", delaylines, dlinenum)
        dlinenum, daystr = getToLine("START DAY", delaylines, dlinenum)
        dlinenum, hourstr = getToLine("START HOUR", delaylines, dlinenum)
        dlinenum, minstr = getToLine("START MINUTE", delaylines, dlinenum)
        dlinenum, secstr = getToLine("START SECOND", delaylines, dlinenum)
        self.startmjd = ymd2mjd(int(yearstr), int(monthstr), int(daystr))
        self.startsec = int(hourstr)*3600 + int(minstr)*60 + int(secstr)
        self.polyorder = 2
        self.polyinterval = 1
        ulinenum, numtelstr = getToLine("NUM TELESCOPES", uvwlines, ulinenum)
        self.numantennas = int(numtelstr)
        self.antennas = []
        self.antennamap = {}
        for i in range(self.numantennas):
            ulinenum, val = getToLine("TELESCOPE %d NAME" % i, uvwlines, ulinenum)
            antennaname = val.strip()
            ulinenum, val = getToLine("TELESCOPE %d X (m)" % i, uvwlines, ulinenum)
            x = float(val.strip())
            ulinenum, val = getToLine("TELESCOPE %d Y (m)" % i, uvwlines, ulinenum)
            y = float(val.strip())
            ulinenum, val = getToLine("TELESCOPE %d Z (m)" % i, uvwlines, ulinenum)
            z = float(val.strip())
            self.antennas.append(Telescope(antennaname, x, y, z, 0.0)) #Lucky the axis offset isn't needed
            self.antennamap[antennaname] = i
        self.sources = []
        dlinenum, val = getToLine("NUM SCANS", delaylines, dlinenum)
        self.numscans = int(val)
        ulinenum, val = getToLine("NUM SCANS", uvwlines, ulinenum)
        if not self.numscans == int(val):
            print "Error - .delay and .uvw file disagree on number of scans!"
            sys.exit(1)
        dlinenum += 1
        ulinenum += 1
        self.scans = []
        for j in range(self.numscans):
            toadd = Scan(delaylines, uvwlines, dlinenum, ulinenum,
                         self.numantennas, self.polyorder, self.polyinterval,
                         self.isDifx15, self.sources, self.startmjd, self.startsec,
                         self.compatibilityMode)
            self.scans.append(toadd)
            dlinenum = toadd.enddlinenum
            ulinenum = toadd.endulinenum
        self.lastscannum = 0
        self.changed = True

    def getDelayAndRate(self, mjd, second, antennaname, sourcename):
        while self.lastscannum < self.numscans and \
              not self.scans[self.lastscannum].containsTime(mjd, second):
            self.lastscannum += 1
            self.changed = True
        if self.lastscannum == self.numscans:
            self.lastscannum = 0
            while self.lastscannum < self.numscans and \
                  not self.scans[self.lastscannum].containsTime(mjd, second):
                self.lastscannum += 1
                self.changed = True
            if self.lastscannum == self.numscans:
                print "Couldn't find scan for time %d/%f" % (mjd, second)
                return -9e99,-9e99
        if self.changed:
            #if self.compatibilityMode:
            #    print "Running in compatibility mode"
            #print "Looking for a source at time " + str(mjd) + "/" + str(second)
            self.sourceindex = self.scans[self.lastscannum].getSourceIndex(sourcename)
            if self.sourceindex < -1: #Couldn't find a source and should have been able to
                print "Was looking for a source at time " + str(mjd) + "/" + str(second)
                print "Current scan runs from " + str(self.scans[self.lastscannum].startmjd) + " + " + str(self.scans[self.lastscannum].scanstart) + \
                      " to " +  str(self.scans[self.lastscannum].scanstart + self.scans[self.lastscannum].scandur)
                self.lastscannum += 1
                if self.scans[self.lastscannum].containsTime(mjd, second):
                    self.sourceindex = self.scans[self.lastscannum].getSourceIndex(sourcename)
                    if self.sourceindex < -1: #Still couldn't find a source and should have been able to
                        print "Still can't find right source - aborting"
                        sys.exit()
                else:
                    print "The next scan ran from " + str(self.scans[self.lastscannum].startmjd) + " + " + str(self.scans[self.lastscannum].scanstart) + \
                      " to " +  str(self.scans[self.lastscannum].scanstart + self.scans[self.lastscannum].scandur)
                    sys.exit()
            self.changed = False
        try:
            antennaindex = self.antennamap[antennaname.strip()]
        except KeyError:
            print "Could not find antenna " + antennaname.strip() + \
                  " in the antenna map from file " + self.filename
            sys.exit()
        if self.sourceindex < 0:
            return -9e99,-9e99
        return self.scans[self.lastscannum].getDelayAndRate(mjd, second, antennaindex,
                                                       self.sourceindex)
    def getUVW(self, antenna1name, antenna2name, sourcename, mjd, second):
        while self.lastscannum < self.numscans and \
              not self.scans[self.lastscannum].containsTime(mjd, second):
            self.lastscannum += 1
            self.changed = True
        if self.lastscannum == self.numscans:
            self.lastscannum = 0
            while self.lastscannum < self.numscans and \
                  not self.scans[self.lastscannum].containsTime(mjd, second):
                self.lastscannum += 1
                self.changed = True
            if self.lastscannum == self.numscans:
                print "Couldn't find scan for time %d/%f" % (mjd, second)
                return [-9e99,-9e99,-9e99]
        if self.changed:
            self.sourceindex = self.scans[self.lastscannum].getSourceIndex(sourcename)
            self.changed = False
        antenna1index = self.antennamap[antenna1name]
        antenna2index = self.antennamap[antenna2name]
        if self.sourceindex < 0:
            return [-9e99,-9e99,-9e99]
        return self.scans[self.lastscannum].getUVW(mjd, second, antenna1index, antenna2index,
                                              self.sourceindex)

