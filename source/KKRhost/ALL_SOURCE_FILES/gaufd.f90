!     ************************************************
SUBROUTINE gaufd(xi,wi,n)
!     ************************************************
!     .. Scalar Arguments ..
INTEGER :: n
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION :: wi(*),xi(*)
!     ..
IF (n == 1) THEN
  xi(1) = -49817229548128141768.d-20
  wi(1) = 10000000000000031192.d-19
  GO TO 10
  
END IF

IF (n == 2) THEN
  xi(1) = -78465071850839016234.d-20
  xi(2) = -20091536266094051757.d-20
  wi(1) = 50923235990870048433.d-20
  wi(2) = 49076764009130263488.d-20
  GO TO 10
  
END IF

IF (n == 3) THEN
  xi(1) = -88288518955458358024.d-20
  xi(2) = -48117621892777473749.d-20
  xi(3) = -88198184413497647625.d-21
  wi(1) = 28858444436509900908.d-20
  wi(2) = 45966895698954759346.d-20
  wi(3) = 25174659864535651667.d-20
  GO TO 10
  
END IF

IF (n == 4) THEN
  xi(1) = -92613063531202843773.d-20
  xi(2) = -64918327008663578157.d-20
  xi(3) = -28982568853420020298.d-20
  xi(4) = -24595209663255169680.d-21
  wi(1) = 18501429405165520392.d-20
  wi(2) = 34614391006511784214.d-20
  wi(3) = 34152482191988153127.d-20
  wi(4) = 12731697396334854188.d-20
  GO TO 10
  
END IF

IF (n == 5) THEN
  xi(1) = -94875333872503463082.d-20
  xi(2) = -74805843506753178608.d-20
  xi(3) = -45504655263391074765.d-20
  xi(4) = -16657582360358973599.d-20
  xi(5) = 27402283545708211900.d-21
  wi(1) = 12939804504572789754.d-20
  wi(2) = 26102400189213290231.d-20
  wi(3) = 30851911091450589451.d-20
  wi(4) = 24746815229701880449.d-20
  wi(5) = 53590689850617620359.d-21
  GO TO 10
  
END IF

IF (n == 6) THEN
  xi(1) = -96204950250095729781.d-20
  xi(2) = -80971428101130972258.d-20
  xi(3) = -57293627456482418171.d-20
  xi(4) = -30755197635518367504.d-20
  xi(5) = -82123839469384988331.d-21
  xi(6) = 83748358371240941581.d-21
  wi(1) = 96268650841705383829.d-21
  wi(2) = 20246201047059595265.d-20
  wi(3) = 26160719441051813381.d-20
  wi(4) = 25781980698475975536.d-20
  wi(5) = 16683001513553609336.d-20
  wi(6) = 15012322156887800205.d-21
  GO TO 10
  
END IF

IF (n == 7) THEN
  xi(1) = -97053934379083423143.d-20
  xi(2) = -85045695849615413757.d-20
  xi(3) = -65665104053460540522.d-20
  xi(4) = -42357896269371657364.d-20
  xi(5) = -19472732441816555564.d-20
  xi(6) = -19669621223691542539.d-21
  xi(7) = 15142830586888806919.d-20
  wi(1) = 74948008822570509041.d-21
  wi(2) = 16170863905729061704.d-20
  wi(3) = 22007120289205973485.d-20
  wi(4) = 23880411919774885893.d-20
  wi(5) = 20952460047488907594.d-20
  wi(6) = 92465405554445737538.d-21
  wi(7) = 24780240009985858690.d-22
  GO TO 10
  
END IF

IF (n == 8) THEN
  xi(1) = -97630544447925725992.d-20
  xi(2) = -87873822716479965943.d-20
  xi(3) = -71736329217593360204.d-20
  xi(4) = -51463306578144813387.d-20
  xi(5) = -29967081434747298359.d-20
  xi(6) = -10763455942936048359.d-20
  xi(7) = 35963113675701677498.d-21
  xi(8) = 23003149140664609750.d-20
  wi(1) = 60394634019629989770.d-21
  wi(2) = 13252509350880929004.d-20
  wi(3) = 18643612522057003210.d-20
  wi(4) = 21413715867867937533.d-20
  wi(5) = 21005092708864293339.d-20
  wi(6) = 16003068683842947897.d-20
  wi(7) = 36159126989806650464.d-21
  wi(8) = 26624765543536915040.d-23
  GO TO 10
  
END IF

IF (n == 9) THEN
  xi(1) = -98041275487012188695.d-20
  xi(2) = -89918326179154863440.d-20
  xi(3) = -76254129548477842110.d-20
  xi(4) = -58579104527384144901.d-20
  xi(5) = -38924212142470946276.d-20
  xi(6) = -19724340764961096691.d-20
  xi(7) = -40039281758884590381.d-21
  xi(8) = 97228170103579374416.d-21
  xi(9) = 31678885353558278864.d-20
  wi(1) = 49992516372028853833.d-21
  wi(2) = 11099301824870447793.d-20
  wi(3) = 15971411690431220541.d-20
  wi(4) = 19037877203046567198.d-20
  wi(5) = 19869087157813151863.d-20
  wi(6) = 17972334325952047726.d-20
  wi(7) = 10203571121909080322.d-20
  wi(8) = 84501828581921130722.d-22
  wi(9) = 21467529556997868476.d-24
  GO TO 10
  
END IF

IF (n == 10) THEN
  xi(1) = -98345122025502045873.d-20
  xi(2) = -91446749996879318119.d-20
  xi(3) = -79700500547314513626.d-20
  xi(4) = -64189534981349313375.d-20
  xi(5) = -46376588343242516012.d-20
  xi(6) = -28030431525349494354.d-20
  xi(7) = -11327091328726333942.d-20
  xi(8) = 17437648086722052805.d-21
  xi(9) = 16877498338102917782.d-20
  xi(10) = 40960465258252015313.d-20
  wi(1) = 42278597323639457484.d-21
  wi(2) = 94666349251635366832.d-21
  wi(3) = 13843777024241956101.d-20
  wi(4) = 16932936699837666261.d-20
  wi(5) = 18398357022114735352.d-20
  wi(6) = 17939886390638648260.d-20
  wi(7) = 14468854182396060463.d-20
  wi(8) = 46026485095922891703.d-21
  wi(9) = 11890402956686871419.d-22
  wi(10) = 14148408460516817666.d-25
  GO TO 10
  
END IF

IF (n == 11) THEN
  xi(1) = -98576901837451635280.d-20
  xi(2) = -92621727156102677473.d-20
  xi(3) = -82389243156123939088.d-20
  xi(4) = -68670708816882492198.d-20
  xi(5) = -52549052940365991088.d-20
  xi(6) = -35349156561982307316.d-20
  xi(7) = -18652071146560858606.d-20
  xi(8) = -45389164233559550280.d-21
  xi(9) = 76984180593432347734.d-21
  xi(10) = 24899533750455431614.d-20
  xi(11) = 50711636785486806957.d-20
  wi(1) = 36383684790132198923.d-21
  wi(2) = 81985364434128201418.d-21
  wi(3) = 12133566247788805356.d-20
  wi(4) = 15122112006362489825.d-20
  wi(5) = 16900090791849557413.d-20
  wi(6) = 17240157268363508589.d-20
  wi(7) = 15745585899461757802.d-20
  wi(8) = 97600157144810676257.d-21
  wi(9) = 12496828256639735424.d-21
  wi(10) = 11876318920871395759.d-23
  wi(11) = 80046822403386311030.d-27
  GO TO 10
  
END IF

IF (n == 12) THEN
  xi(1) = -98758247347129831371.d-20
  xi(2) = -93546465146779806654.d-20
  xi(3) = -84528996754470930223.d-20
  xi(4) = -72299594230844519839.d-20
  xi(5) = -57679398168141327066.d-20
  xi(6) = -41683730779892996801.d-20
  xi(7) = -25514627335790291149.d-20
  xi(8) = -10710838211747769681.d-20
  xi(9) = 12720145729326415607.d-21
  xi(10) = 14540842218988328389.d-20
  xi(11) = 33552500235752414908.d-20
  xi(12) = 60838109964484063119.d-20
  wi(1) = 31765161579790701148.d-21
  wi(2) = 71927618746964313778.d-21
  wi(3) = 10742555378156694842.d-20
  wi(4) = 13578811351554214795.d-20
  wi(5) = 15492042553417744038.d-20
  wi(6) = 16300300254834219520.d-20
  wi(7) = 15784577013790806216.d-20
  wi(8) = 12921482926208917372.d-20
  wi(9) = 46096943233133302568.d-21
  wi(10) = 20030610755774790850.d-22
  wi(11) = 95165705752725893549.d-25
  wi(12) = 40143360822128708729.d-28
  GO TO 10
  
END IF

IF (n == 13) THEN
  xi(1) = -98903182721370020265.d-20
  xi(2) = -94288936524363459773.d-20
  xi(3) = -86261843870640242196.d-20
  xi(4) = -75277808759167753869.d-20
  xi(5) = -61972590294795871779.d-20
  xi(6) = -47139332563986024748.d-20
  xi(7) = -31718188942187627557.d-20
  xi(8) = -16854863011308355787.d-20
  xi(9) = -41195843159851553906.d-21
  xi(10) = 71957380142115164738.d-21
  xi(11) = 22223926926874000328.d-20
  xi(12) = 42682885634093164862.d-20
  xi(13) = 71270930856714354732.d-20
  wi(1) = 28069991026027589482.d-21
  wi(2) = 63803895087070663653.d-21
  wi(3) = 95973484361405430270.d-21
  wi(4) = 12264378189747678145.d-20
  wi(5) = 14213612346123977130.d-20
  wi(6) = 15296686007570952707.d-20
  wi(7) = 15358437552921000921.d-20
  wi(8) = 14007635729175637795.d-20
  wi(9) = 87531230524252970103.d-21
  wi(10) = 12989730151883234012.d-21
  wi(11) = 22351943999969127535.d-23
  wi(12) = 65097139765619073344.d-26
  wi(13) = 18257341724040876662.d-29
  GO TO 10
  
END IF

IF (n == 14) THEN
  xi(1) = -99021130855943209687.d-20
  xi(2) = -94895368426058288869.d-20
  xi(3) = -87686856465753704289.d-20
  xi(4) = -77752669471002194917.d-20
  xi(5) = -65594116901081876554.d-20
  xi(6) = -51841232227159879604.d-20
  xi(7) = -37243750660439082187.d-20
  xi(8) = -22693429290756856295.d-20
  xi(9) = -93940943648510570987.d-21
  xi(10) = 16521198218716065629.d-21
  xi(11) = 13919799114797561344.d-20
  xi(12) = 30521886852802066309.d-20
  xi(13) = 52192337126752562221.d-20
  xi(14) = 81957965081548293179.d-20
  wi(1) = 25060310888021301605.d-21
  wi(2) = 57137272611562033779.d-21
  wi(3) = 86434450014324433897.d-21
  wi(4) = 11141118228632175288.d-20
  wi(5) = 13070790263291078499.d-20
  wi(6) = 14310195071194851995.d-20
  wi(7) = 14737968606274298328.d-20
  wi(8) = 14154903694980505066.d-20
  wi(9) = 11456160782223814050.d-20
  wi(10) = 40466499493397342820.d-21
  wi(11) = 21701008894932486895.d-22
  wi(12) = 19960253076851250807.d-24
  wi(13) = 39376501060604877095.d-27
  wi(14) = 76596142918862399780.d-31
  GO TO 10
  
END IF

IF (n == 15) THEN
  xi(1) = -99118619138431485634.d-20
  xi(2) = -95398089203095832045.d-20
  xi(3) = -88874665207045485764.d-20
  xi(4) = -79832886799647722652.d-20
  xi(5) = -68674462209286747178.d-20
  xi(6) = -55907326778454372362.d-20
  xi(7) = -42138595122137487519.d-20
  xi(8) = -28083407355763995168.d-20
  xi(9) = -14649293944496725019.d-20
  xi(10) = -30865949117072113052.d-21
  xi(11) = 75989566859912966734.d-21
  xi(12) = 21425891814116860148.d-20
  xi(13) = 39280262275215780450.d-20
  xi(14) = 62012182191671475949.d-20
  xi(15) = 92858877219218103945.d-20
  wi(1) = 22570991165870390473.d-21
  wi(2) = 51589746641923392000.d-21
  wi(3) = 78401918844466166239.d-21
  wi(4) = 10176234626640128024.d-20
  wi(5) = 12055819130110177262.d-20
  wi(6) = 13377324647273569326.d-20
  wi(7) = 14041818603247829422.d-20
  wi(8) = 13919569003129657925.d-20
  wi(9) = 12562361445602688222.d-20
  wi(10) = 74852662340708470150.d-21
  wi(11) = 10996744175647251144.d-21
  wi(12) = 25513307315040157893.d-23
  wi(13) = 15270418102934789627.d-25
  wi(14) = 21560859319293022163.d-28
  wi(15) = 30032040385910287756.d-32
  GO TO 10
  
END IF

IF (n == 16) THEN
  xi(1) = -99200289748411473927.d-20
  xi(2) = -95820266378296634182.d-20
  xi(3) = -89876661129475763142.d-20
  xi(4) = -81599671254865832971.d-20
  xi(5) = -71315812647978444249.d-20
  xi(6) = -59440032425488487666.d-20
  xi(7) = -46470396871945791541.d-20
  xi(8) = -32991653294098863600.d-20
  xi(9) = -19716091326862980561.d-20
  xi(10) = -76605243443508959615.d-21
  xi(11) = 26155046503992069925.d-21
  xi(12) = 14307776307824938682.d-20
  xi(13) = 29506185654032182160.d-20
  xi(14) = 48403577800553841578.d-20
  xi(15) = 72091584865612160132.d-20
  xi(16) = 10394188783181811718.d-19
  wi(1) = 20484388078614008045.d-21
  wi(2) = 46916532350372347409.d-21
  wi(3) = 71569877291069983495.d-21
  wi(4) = 93424466379672137196.d-21
  wi(5) = 11156011364306951297.d-20
  wi(6) = 12512553084306063601.d-20
  wi(7) = 13329704953113185969.d-20
  wi(8) = 13510959073859290681.d-20
  wi(9) = 12840858805365359846.d-20
  wi(10) = 10016528657871746742.d-20
  wi(11) = 32102655847303900301.d-21
  wi(12) = 18115418480524121495.d-22
  wi(13) = 24274994772381143993.d-24
  wi(14) = 10371321943363515335.d-26
  wi(15) = 10868941709467004901.d-29
  wi(16) = 11117372791599461059.d-33
  GO TO 10
  
END IF

10 CONTINUE
RETURN

END !       SUBROUTINE GAUFD(XI,WI,N)
