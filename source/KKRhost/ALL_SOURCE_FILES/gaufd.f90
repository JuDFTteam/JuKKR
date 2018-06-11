!     ************************************************
    Subroutine gaufd(xi, wi, n)
      Use mod_datatypes, Only: dp
!     ************************************************
!     .. Scalar Arguments ..
      Integer :: n
!     ..
!     .. Array Arguments ..
      Real (Kind=dp) :: wi(*), xi(*)
!     ..
      If (n==1) Then
        xi(1) = -49817229548128141768.E-20_dp
        wi(1) = 10000000000000031192.E-19_dp
        Go To 100

      End If

      If (n==2) Then
        xi(1) = -78465071850839016234.E-20_dp
        xi(2) = -20091536266094051757.E-20_dp
        wi(1) = 50923235990870048433.E-20_dp
        wi(2) = 49076764009130263488.E-20_dp
        Go To 100

      End If

      If (n==3) Then
        xi(1) = -88288518955458358024.E-20_dp
        xi(2) = -48117621892777473749.E-20_dp
        xi(3) = -88198184413497647625.E-21_dp
        wi(1) = 28858444436509900908.E-20_dp
        wi(2) = 45966895698954759346.E-20_dp
        wi(3) = 25174659864535651667.E-20_dp
        Go To 100

      End If

      If (n==4) Then
        xi(1) = -92613063531202843773.E-20_dp
        xi(2) = -64918327008663578157.E-20_dp
        xi(3) = -28982568853420020298.E-20_dp
        xi(4) = -24595209663255169680.E-21_dp
        wi(1) = 18501429405165520392.E-20_dp
        wi(2) = 34614391006511784214.E-20_dp
        wi(3) = 34152482191988153127.E-20_dp
        wi(4) = 12731697396334854188.E-20_dp
        Go To 100

      End If

      If (n==5) Then
        xi(1) = -94875333872503463082.E-20_dp
        xi(2) = -74805843506753178608.E-20_dp
        xi(3) = -45504655263391074765.E-20_dp
        xi(4) = -16657582360358973599.E-20_dp
        xi(5) = 27402283545708211900.E-21_dp
        wi(1) = 12939804504572789754.E-20_dp
        wi(2) = 26102400189213290231.E-20_dp
        wi(3) = 30851911091450589451.E-20_dp
        wi(4) = 24746815229701880449.E-20_dp
        wi(5) = 53590689850617620359.E-21_dp
        Go To 100

      End If

      If (n==6) Then
        xi(1) = -96204950250095729781.E-20_dp
        xi(2) = -80971428101130972258.E-20_dp
        xi(3) = -57293627456482418171.E-20_dp
        xi(4) = -30755197635518367504.E-20_dp
        xi(5) = -82123839469384988331.E-21_dp
        xi(6) = 83748358371240941581.E-21_dp
        wi(1) = 96268650841705383829.E-21_dp
        wi(2) = 20246201047059595265.E-20_dp
        wi(3) = 26160719441051813381.E-20_dp
        wi(4) = 25781980698475975536.E-20_dp
        wi(5) = 16683001513553609336.E-20_dp
        wi(6) = 15012322156887800205.E-21_dp
        Go To 100

      End If

      If (n==7) Then
        xi(1) = -97053934379083423143.E-20_dp
        xi(2) = -85045695849615413757.E-20_dp
        xi(3) = -65665104053460540522.E-20_dp
        xi(4) = -42357896269371657364.E-20_dp
        xi(5) = -19472732441816555564.E-20_dp
        xi(6) = -19669621223691542539.E-21_dp
        xi(7) = 15142830586888806919.E-20_dp
        wi(1) = 74948008822570509041.E-21_dp
        wi(2) = 16170863905729061704.E-20_dp
        wi(3) = 22007120289205973485.E-20_dp
        wi(4) = 23880411919774885893.E-20_dp
        wi(5) = 20952460047488907594.E-20_dp
        wi(6) = 92465405554445737538.E-21_dp
        wi(7) = 24780240009985858690.E-22_dp
        Go To 100

      End If

      If (n==8) Then
        xi(1) = -97630544447925725992.E-20_dp
        xi(2) = -87873822716479965943.E-20_dp
        xi(3) = -71736329217593360204.E-20_dp
        xi(4) = -51463306578144813387.E-20_dp
        xi(5) = -29967081434747298359.E-20_dp
        xi(6) = -10763455942936048359.E-20_dp
        xi(7) = 35963113675701677498.E-21_dp
        xi(8) = 23003149140664609750.E-20_dp
        wi(1) = 60394634019629989770.E-21_dp
        wi(2) = 13252509350880929004.E-20_dp
        wi(3) = 18643612522057003210.E-20_dp
        wi(4) = 21413715867867937533.E-20_dp
        wi(5) = 21005092708864293339.E-20_dp
        wi(6) = 16003068683842947897.E-20_dp
        wi(7) = 36159126989806650464.E-21_dp
        wi(8) = 26624765543536915040.E-23_dp
        Go To 100

      End If

      If (n==9) Then
        xi(1) = -98041275487012188695.E-20_dp
        xi(2) = -89918326179154863440.E-20_dp
        xi(3) = -76254129548477842110.E-20_dp
        xi(4) = -58579104527384144901.E-20_dp
        xi(5) = -38924212142470946276.E-20_dp
        xi(6) = -19724340764961096691.E-20_dp
        xi(7) = -40039281758884590381.E-21_dp
        xi(8) = 97228170103579374416.E-21_dp
        xi(9) = 31678885353558278864.E-20_dp
        wi(1) = 49992516372028853833.E-21_dp
        wi(2) = 11099301824870447793.E-20_dp
        wi(3) = 15971411690431220541.E-20_dp
        wi(4) = 19037877203046567198.E-20_dp
        wi(5) = 19869087157813151863.E-20_dp
        wi(6) = 17972334325952047726.E-20_dp
        wi(7) = 10203571121909080322.E-20_dp
        wi(8) = 84501828581921130722.E-22_dp
        wi(9) = 21467529556997868476.E-24_dp
        Go To 100

      End If

      If (n==10) Then
        xi(1) = -98345122025502045873.E-20_dp
        xi(2) = -91446749996879318119.E-20_dp
        xi(3) = -79700500547314513626.E-20_dp
        xi(4) = -64189534981349313375.E-20_dp
        xi(5) = -46376588343242516012.E-20_dp
        xi(6) = -28030431525349494354.E-20_dp
        xi(7) = -11327091328726333942.E-20_dp
        xi(8) = 17437648086722052805.E-21_dp
        xi(9) = 16877498338102917782.E-20_dp
        xi(10) = 40960465258252015313.E-20_dp
        wi(1) = 42278597323639457484.E-21_dp
        wi(2) = 94666349251635366832.E-21_dp
        wi(3) = 13843777024241956101.E-20_dp
        wi(4) = 16932936699837666261.E-20_dp
        wi(5) = 18398357022114735352.E-20_dp
        wi(6) = 17939886390638648260.E-20_dp
        wi(7) = 14468854182396060463.E-20_dp
        wi(8) = 46026485095922891703.E-21_dp
        wi(9) = 11890402956686871419.E-22_dp
        wi(10) = 14148408460516817666.E-25_dp
        Go To 100

      End If

      If (n==11) Then
        xi(1) = -98576901837451635280.E-20_dp
        xi(2) = -92621727156102677473.E-20_dp
        xi(3) = -82389243156123939088.E-20_dp
        xi(4) = -68670708816882492198.E-20_dp
        xi(5) = -52549052940365991088.E-20_dp
        xi(6) = -35349156561982307316.E-20_dp
        xi(7) = -18652071146560858606.E-20_dp
        xi(8) = -45389164233559550280.E-21_dp
        xi(9) = 76984180593432347734.E-21_dp
        xi(10) = 24899533750455431614.E-20_dp
        xi(11) = 50711636785486806957.E-20_dp
        wi(1) = 36383684790132198923.E-21_dp
        wi(2) = 81985364434128201418.E-21_dp
        wi(3) = 12133566247788805356.E-20_dp
        wi(4) = 15122112006362489825.E-20_dp
        wi(5) = 16900090791849557413.E-20_dp
        wi(6) = 17240157268363508589.E-20_dp
        wi(7) = 15745585899461757802.E-20_dp
        wi(8) = 97600157144810676257.E-21_dp
        wi(9) = 12496828256639735424.E-21_dp
        wi(10) = 11876318920871395759.E-23_dp
        wi(11) = 80046822403386311030.E-27_dp
        Go To 100

      End If

      If (n==12) Then
        xi(1) = -98758247347129831371.E-20_dp
        xi(2) = -93546465146779806654.E-20_dp
        xi(3) = -84528996754470930223.E-20_dp
        xi(4) = -72299594230844519839.E-20_dp
        xi(5) = -57679398168141327066.E-20_dp
        xi(6) = -41683730779892996801.E-20_dp
        xi(7) = -25514627335790291149.E-20_dp
        xi(8) = -10710838211747769681.E-20_dp
        xi(9) = 12720145729326415607.E-21_dp
        xi(10) = 14540842218988328389.E-20_dp
        xi(11) = 33552500235752414908.E-20_dp
        xi(12) = 60838109964484063119.E-20_dp
        wi(1) = 31765161579790701148.E-21_dp
        wi(2) = 71927618746964313778.E-21_dp
        wi(3) = 10742555378156694842.E-20_dp
        wi(4) = 13578811351554214795.E-20_dp
        wi(5) = 15492042553417744038.E-20_dp
        wi(6) = 16300300254834219520.E-20_dp
        wi(7) = 15784577013790806216.E-20_dp
        wi(8) = 12921482926208917372.E-20_dp
        wi(9) = 46096943233133302568.E-21_dp
        wi(10) = 20030610755774790850.E-22_dp
        wi(11) = 95165705752725893549.E-25_dp
        wi(12) = 40143360822128708729.E-28_dp
        Go To 100

      End If

      If (n==13) Then
        xi(1) = -98903182721370020265.E-20_dp
        xi(2) = -94288936524363459773.E-20_dp
        xi(3) = -86261843870640242196.E-20_dp
        xi(4) = -75277808759167753869.E-20_dp
        xi(5) = -61972590294795871779.E-20_dp
        xi(6) = -47139332563986024748.E-20_dp
        xi(7) = -31718188942187627557.E-20_dp
        xi(8) = -16854863011308355787.E-20_dp
        xi(9) = -41195843159851553906.E-21_dp
        xi(10) = 71957380142115164738.E-21_dp
        xi(11) = 22223926926874000328.E-20_dp
        xi(12) = 42682885634093164862.E-20_dp
        xi(13) = 71270930856714354732.E-20_dp
        wi(1) = 28069991026027589482.E-21_dp
        wi(2) = 63803895087070663653.E-21_dp
        wi(3) = 95973484361405430270.E-21_dp
        wi(4) = 12264378189747678145.E-20_dp
        wi(5) = 14213612346123977130.E-20_dp
        wi(6) = 15296686007570952707.E-20_dp
        wi(7) = 15358437552921000921.E-20_dp
        wi(8) = 14007635729175637795.E-20_dp
        wi(9) = 87531230524252970103.E-21_dp
        wi(10) = 12989730151883234012.E-21_dp
        wi(11) = 22351943999969127535.E-23_dp
        wi(12) = 65097139765619073344.E-26_dp
        wi(13) = 18257341724040876662.E-29_dp
        Go To 100

      End If

      If (n==14) Then
        xi(1) = -99021130855943209687.E-20_dp
        xi(2) = -94895368426058288869.E-20_dp
        xi(3) = -87686856465753704289.E-20_dp
        xi(4) = -77752669471002194917.E-20_dp
        xi(5) = -65594116901081876554.E-20_dp
        xi(6) = -51841232227159879604.E-20_dp
        xi(7) = -37243750660439082187.E-20_dp
        xi(8) = -22693429290756856295.E-20_dp
        xi(9) = -93940943648510570987.E-21_dp
        xi(10) = 16521198218716065629.E-21_dp
        xi(11) = 13919799114797561344.E-20_dp
        xi(12) = 30521886852802066309.E-20_dp
        xi(13) = 52192337126752562221.E-20_dp
        xi(14) = 81957965081548293179.E-20_dp
        wi(1) = 25060310888021301605.E-21_dp
        wi(2) = 57137272611562033779.E-21_dp
        wi(3) = 86434450014324433897.E-21_dp
        wi(4) = 11141118228632175288.E-20_dp
        wi(5) = 13070790263291078499.E-20_dp
        wi(6) = 14310195071194851995.E-20_dp
        wi(7) = 14737968606274298328.E-20_dp
        wi(8) = 14154903694980505066.E-20_dp
        wi(9) = 11456160782223814050.E-20_dp
        wi(10) = 40466499493397342820.E-21_dp
        wi(11) = 21701008894932486895.E-22_dp
        wi(12) = 19960253076851250807.E-24_dp
        wi(13) = 39376501060604877095.E-27_dp
        wi(14) = 76596142918862399780.E-31_dp
        Go To 100

      End If

      If (n==15) Then
        xi(1) = -99118619138431485634.E-20_dp
        xi(2) = -95398089203095832045.E-20_dp
        xi(3) = -88874665207045485764.E-20_dp
        xi(4) = -79832886799647722652.E-20_dp
        xi(5) = -68674462209286747178.E-20_dp
        xi(6) = -55907326778454372362.E-20_dp
        xi(7) = -42138595122137487519.E-20_dp
        xi(8) = -28083407355763995168.E-20_dp
        xi(9) = -14649293944496725019.E-20_dp
        xi(10) = -30865949117072113052.E-21_dp
        xi(11) = 75989566859912966734.E-21_dp
        xi(12) = 21425891814116860148.E-20_dp
        xi(13) = 39280262275215780450.E-20_dp
        xi(14) = 62012182191671475949.E-20_dp
        xi(15) = 92858877219218103945.E-20_dp
        wi(1) = 22570991165870390473.E-21_dp
        wi(2) = 51589746641923392000.E-21_dp
        wi(3) = 78401918844466166239.E-21_dp
        wi(4) = 10176234626640128024.E-20_dp
        wi(5) = 12055819130110177262.E-20_dp
        wi(6) = 13377324647273569326.E-20_dp
        wi(7) = 14041818603247829422.E-20_dp
        wi(8) = 13919569003129657925.E-20_dp
        wi(9) = 12562361445602688222.E-20_dp
        wi(10) = 74852662340708470150.E-21_dp
        wi(11) = 10996744175647251144.E-21_dp
        wi(12) = 25513307315040157893.E-23_dp
        wi(13) = 15270418102934789627.E-25_dp
        wi(14) = 21560859319293022163.E-28_dp
        wi(15) = 30032040385910287756.E-32_dp
        Go To 100

      End If

      If (n==16) Then
        xi(1) = -99200289748411473927.E-20_dp
        xi(2) = -95820266378296634182.E-20_dp
        xi(3) = -89876661129475763142.E-20_dp
        xi(4) = -81599671254865832971.E-20_dp
        xi(5) = -71315812647978444249.E-20_dp
        xi(6) = -59440032425488487666.E-20_dp
        xi(7) = -46470396871945791541.E-20_dp
        xi(8) = -32991653294098863600.E-20_dp
        xi(9) = -19716091326862980561.E-20_dp
        xi(10) = -76605243443508959615.E-21_dp
        xi(11) = 26155046503992069925.E-21_dp
        xi(12) = 14307776307824938682.E-20_dp
        xi(13) = 29506185654032182160.E-20_dp
        xi(14) = 48403577800553841578.E-20_dp
        xi(15) = 72091584865612160132.E-20_dp
        xi(16) = 10394188783181811718.E-19_dp
        wi(1) = 20484388078614008045.E-21_dp
        wi(2) = 46916532350372347409.E-21_dp
        wi(3) = 71569877291069983495.E-21_dp
        wi(4) = 93424466379672137196.E-21_dp
        wi(5) = 11156011364306951297.E-20_dp
        wi(6) = 12512553084306063601.E-20_dp
        wi(7) = 13329704953113185969.E-20_dp
        wi(8) = 13510959073859290681.E-20_dp
        wi(9) = 12840858805365359846.E-20_dp
        wi(10) = 10016528657871746742.E-20_dp
        wi(11) = 32102655847303900301.E-21_dp
        wi(12) = 18115418480524121495.E-22_dp
        wi(13) = 24274994772381143993.E-24_dp
        wi(14) = 10371321943363515335.E-26_dp
        wi(15) = 10868941709467004901.E-29_dp
        wi(16) = 11117372791599461059.E-33_dp
        Go To 100

      End If

100   Continue
      Return

    End Subroutine !       SUBROUTINE GAUFD(XI,WI,N)
