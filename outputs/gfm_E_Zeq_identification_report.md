# GFMVSG1: E ve Zeq'nin fiziksel tanımlanması

Tarih: 1 Ekim 2026. Bu rapor, model topolojisi ve parametrelerinden türetilen **dengeli, pozitif sıralı, 50 Hz quasi-steady fazör indirgemesidir**. EMT eğri uydurması, parametre optimizasyonu veya yeni contraction hesabı yapılmadı. İki mevcut MATLAB analiz dosyası değiştirilmedi.

## 1. Sonuç ve kapsam

En tutarlı seçim, **gerilim döngüsünün düzenlediği kapasitör/filtre düğümü gerilimidir**. Normal, limiter-inaktif durumda ideal iç-döngü takibi ve ölçümlerde birim kazanç yaklaşımıyla:

\[
E:=|\underline V_c|,\qquad
E=1+0.1(Q^\star-Q_g).
\]

Qstar=0 ve analitik olarak istenen Pstar=1 pu için ağ ve Q–V denklemlerini birlikte çözünce:

| Büyüklük | Değer | Açıklama |
|---|---:|---|
| Nominal E | 0.999003208451 pu | 690 V hatlar arası RMS / 2 MVA bazı |
| Fiziksel sonsuz bara, LV'ye indirgenmiş | 0.980392156863 pu | 20 kV / 20.4 kV |
| Thevenin kaynak büyüklüğü Uth | 0.983680194662 pu | Şöntler dahil port eşdeğeri |
| Req | 0.016884633807 pu | Aşağıda tanımlanan **dönüştürücü akımı portu** eşdeğeri |
| Xeq | 0.066199725744 pu | Aynı port ve baz |
| abs(Zeq) | 0.068319064304 pu | Fiziksel empedans bazı 0.23805 ohm |
| angle(Zeq) | 75.6914383263 deg | Direnç ihmal edilmedi |
| Imax | 1.2 pu | Dönüştürücü akım referansı, 2 MVA / 690 V bazı |
| Nominal gamma, fiziksel bara referanslı | 3.84979056312 deg | Trafo oranı/fazı LV'ye taşındıktan sonra |
| Nominal gamma, Thevenin referanslı | 3.90344727653 deg | Akım çemberi formülünde kullanılan açı |
| Limiter eşiği, nominal E sabit | 4.65622027424 deg | Thevenin referanslı |
| Limiter eşiği, cebirsel Q–V korunursa | 4.62826772873 deg | E eşikte 1.00099055907 pu olur |

**Bu tablo tek başına `Xeq` değiştirilerek eski koda konulmamalıdır.** Uth fiziksel Ubus'tan farklıdır; limiter akımı ile VSM'nin güç ölçüm yeri de farklıdır. E'nin Q'ya bağlılığı ve açı referansı açıkça korunmalıdır.

Bilimsel değerlendirme **B: yaklaşık ama savunulabilir**. Lineer 50 Hz pasif ağın port-Thevenin eşdeğeri seçilen varsayımlar altında tamdır. İç kontrolün ideal izleme olarak indirgenmesi, Q–V dinamiğinin cebirsel hale getirilmesi ve tüm geçişlerin iki durumla temsil edilmesi ise yaklaşımdır. Aşağıdaki sayılar bir EMT çalışma noktası ölçümü değildir.

## 2. Kaynaklar ve fiilen kullanılan yapı

Ana model: `model/GFMVSG1.slx`.

Kütüphane kökü: `model/tdmsquared-master/tdmsquared-master/Library`.

İncelenen dosyalar:

- `Distributed Generation Units/Converter/Grid_Forming_Converter/Model_DG_Grid_Forming_Converter.slx`
- Aynı klasörde `Subsys_GFormConv_Control.slx`, `Subsys_GFormConv_OuterCtrl_VSM.slx`, `Subsys_GFormConv_InnerCtrl_cascadedPI.slx`
- `Converter_Shared_Subsystems/Subsys_Conv_Electrical_Part_SPS.slx`, `Subsys_Conv_Ctrl_measurement_transformations.slx` (Converter klasörü altında)
- `Shared Functions/dq_Limiter.slx`, `PI_freeze.slx`, `DDSRF.slx`

SLX arşivlerindeki bağlantılar, aktif variant parametreleri ve mask initialization ifadeleri okundu. MATLAB R2024b'de model salt okunur olarak yüklendi; VSM kütüphanesinin yukarıdaki ağaca çözümlendiği ve P_set Signal Builder verileri doğrulandı. Model kaydedilmedi; EMT simülasyonu çalıştırılmadı.

Kaydedilmiş modelde **P_set=0.2, 0.7, 0.8 MW**, yani 2 MVA bazında **0.1, 0.35, 0.4 pu** seviyeleri var. Zaman düğümleri `[0,5,10,12,25,27,35]` s ve değerler `[0.2,0.2,0.2,0.7,0.7,0.8,0.8]` MW. Dolayısıyla Pstar=1 pu burada kullanıcının analitik nominal noktasıdır; model dosyasının kayıtlı yükleme profili değildir. Q_set=0, V_set=20 kV, f_set=50 Hz.

Modelin kayıtlı VSM tuning'i H=0.1, K=20, Kd=40; analitik çalışmanın H=3.4, Kd=5 tuning'inden farklıdır. Bu farklılık elektriksel baz/empedans hesabına eklenmedi ve hiçbir değer değiştirilmedi.

## 3. E adaylarının karşılaştırılması

| Aday/sinyal | VSG açısı ve genlik ilişkisi | Konumu ve geride kalan ağ | Karar |
|---|---|---|---|
| VSM `v_grid*`, Q–V çıkışı | Genlik komutu; normal modda VSM theta ile birlikte kullanılır | İç döngülerden önce; fiziksel bir port değildir | E'nin referansını verir, tek başına fiziksel E değildir |
| `Operation_Mode_Selector/v_grid*_d,q` | Lvirt=0 ile d=Vouter, q=0; komut açısı theta | Kapasitör gerilim hedefi | İdeal takip varsayımının doğru başlangıcı |
| Fiziksel `V_cap_abc` | Limiter-inaktif ideal izleme altında referansla eşleşir; ölçüm fazı düzeltilebilir | Dönüştürücü RL sonrasında, trafo öncesinde | **Seçilen fiziksel E düğümü** |
| Ölçülen `dq+.v_cap_dq+` | PI'nin kontrol ettiği filtrelenmiş gerilim; sıfır hatada `[Eref,0]` | Fiziksel düğümün ölçüm temsili | Fiziksel gerilim ile filtre transferi üzerinden ilişkilendirilmeli |
| Akım PI ham çıkışı `u_d,u_q` | İki bileşen de akım, integratör ve decoupling'e bağlı; arg(u) genel olarak sıfır değil | 1.3 pu limiter öncesinde | Açısı yalnız theta değildir; E olarak seçilmez |
| Gerilim limiteri çıkışı `w` | Radyal sınırlandırılmış u; açı offset'i devam eder | Aktif sönümleme öncesinde | 1.3 pu bir nominal E değeri değildir |
| Aktif sönümleme sonrası `v_cmd` | `w-0.2 i_conv,dq0`; akıma bağlı vektör | Köprü giriş komutu | Sabit genlikli VSG kaynağı değildir |
| Fiziksel köprü `V_conv_abc` | Komut ve dönüştürücüye bağlı; fiziksel açısı theta'ya ek yük açısı içerir | Önünde converter RL ve şönt RC de kalır | Alternatif elektrik portu mümkündür, fakat aynı gamma ile basit VSG kaynağı olmaz |
| Trafo sonrası `V_grid_abc` | Ağ tarafından belirlenir; doğrudan d=E,q=0 kontrol hedefi değil | Hat ve kaynak empedansı geride kalır | Güç ölçüm yeri; iç kaynak E'si değildir |

Sinyal zinciri: `Outer_Control/Virtual_Synchronous_Machine` Q–V çıkışı → `Operation_Mode_Selector/Virtual_Inductance` (Lvirt=0) → `Inner_Control/Cascaded_PI_controller/v_cap_Controller` → akım limiteri → akım PI → gerilim limiteri → Active_Damping → ters Park → köprü → converter RL → **V_cap düğümü**.

Gerilim PI integratörleri etkin ve sınırlayıcılar inaktifken kalıcı dq hata sıfır olmalıdır. Bu nedenle seçilen düğüm için nominal takip yaklaşımı fiziksel olarak anlamlıdır. Freeze/limiter aktifken bu koşul genel olarak geçmez; sonuçlarımız ilk aktivasyonun normal çalışma tarafına ilişkindir.

### E neden yaklaşık 1 pu?

`setpoints_SI_to_pu/Gain2`, V_set değerini **20 kV nominal grid gerilimine** böler: 20/20=1. VSM'de Q–V DC kazancı 0.1; TQ=0.2 s filtresinin DC kazancı 1'dir. Normal modda Lvirt=0 olduğu için:

\[
E_{ref}=1+0.1(0-Q_{meas}).
\]

Birim ölçüm/takip yaklaşımıyla P_g=1 ve ağ denklemleri birlikte çözülünce Q_g=0.009967915484 pu; buradan E=0.999003208451 pu çıkar. **Qstar=0 olması, Q_g'nin sıfır olduğu anlamına gelmez.** Q–V droop integral Q regülatörü değildir.

Bu değeri yerel nominal analizde sabit almak makuldür: nominal noktadan ilk akım limitine kadar E değişimi yaklaşık 0.199%'dur. Tüm açılar, sag'ler veya limiter-aktif çalışma için sabitlik kanıtlanmış değildir. İki durum korunarak E'yi cebirsel Q–V denklemiyle çözmek daha tutarlı seçenektir.

## 4. Gerçek elektrik topolojisi

Pozitif sıra, normal kapalı kesici için:

```text
bridge -- Zf -- Vc -- Zt1 -- core node -- ideal transformer / Zt2 -- HV node
                 |              |                                  |
              series RC        Ym                               breaker
                 |              |                                  |
               ground         ground                     grid power measurement
                                                                  |
                                                         line pi -- Zsource -- Uinf
```

Trafo oranı ve fazı LV'ye taşınınca ideal transformer diyagramdan kaldırılabilir. Trafo **Yg / Delta(D1)**, nominal oranı **690 / (1.02*20000)=690/20400**. D1 gerilimi yıldız tarafına göre 30° geridedir; dolayısıyla HV sonsuz bara fazörü LV'ye taşınırken +30° dönüş uygulanır. Bu bağlantı anlamı [MathWorks transformer dokümanında](https://jp.mathworks.com/help/simscape-electrical/powersys/ref/threephasetransformertwowindings.html) da belirtilir.

Fiziksel gamma, `arg(Vc)-arg(Uinf_referred)` olarak tanımlanır. Doğrudan `theta_VSM-theta_grid_HV` aynı açı değildir. Thevenin gamma'sı ise buna ek olarak aşağıdaki Uth'nin küçük fazını çıkarır.

### Dahil edilen bileşenler

Tüm pu değerleri 2 MVA ve ideal trafoyla uyumlu 690/20400 V bazındadır:

| Eleman | Kullanılan değer | Gerekçe |
|---|---:|---|
| Trafo LV kaçağı Zt1 | 0.005+j0.02 pu | Vc ile core arası |
| Trafo HV kaçağı Zt2 | 0.005+j0.02 pu | HV winding kendi 20.4 kV nominal bazında |
| Mıknatıslanma Ym | 0.002−j0.002 pu | Rm=Lm=500; doygunluk kapalı |
| Kapasitör+damping kolu Yc | 0.000143519441+j0.045689549179 pu | `1/(0.06875−j/0.04569)` |
| Kapalı kesici | 0.01 ohm || 1 Mohm; 0.000048058439 pu | Aktif HV seri yol |
| Hat seri Zline | 0.004815455594+j0.006881670182 pu | 2 km: R=1.002 ohm, L=4.558 mH |
| Hat her uç şöntü | j0.009883983302 pu | Toplam 302.4 nF'nin yarısı |
| Kaynak Zsource | 0.001912797367+j0.019127973668 pu | 100 MVA kısa-devre gücü ve X/R=10 |
| Açık fault snubber | 0.00020808 pu iletkenlik | Pozitif sırada faz başına 1 Mohm, normal durumda |

Hat PI topolojisinin seri ve şönt ayrımı blok dosyasından alındı. Genel PI hat yapısı için ayrıca [MathWorks hat modeli açıklamasına](https://www.mathworks.com/help/simscape-electrical/ref/transmissionlinethreephase.html) bakılabilir. Şöntler seri reaktansa eklenmedi.

Kaynakta `SpecifyImpedance=on` mask seçeneği kısa-devre seviyesinden türetmeyi etkinleştiriyor. |Zs|=20kV²/100MVA=4 ohm, Rs=4/sqrt(101)=0.398014876084 ohm, Xs=3.980148760840 ohm. Kayıtlı alternatif 0.8929 ohm/16.58 mH alanları bu ayarın aktif hesabı değildir.

### Dışarıda bırakılanlar

- Converter RL `0.001+j0.05097 pu`, **seçilen Vc düğümünün gerisinde** kaldığından Vc-to-grid Zeq'ye eklenmedi. Köprü gerilimi ve limiter başlığı kontrolünde tekrar kullanıldı.
- `L_filter_grid_pu=[0.01,0.04]`, **No_Transformer** alternatifinde. Aktif trafo konfigürasyonuna eklenmedi.
- Aktif sönümleme 0.2 pu fiziksel seri direnç değildir. PI, gerekli gerilim komutunu değiştirerek normal rejimde fiziksel port gerilimini koruyabilir; 0.2'yi pasif ağa eklemek kayıpları ve yük akışını yanlış değiştirir.
- OuterCtrl_Lvirt=0; bu yoldan eklenecek sanal endüktans yok.
- Akım/gerilim limiterleri, PI/freeze durumları ve Zvi/Zlim pasif Zeq'ye katılmadı.
- Zero-sequence/core üç-bacak ayrıntıları dengeli pozitif sıra çalışmasında kullanılmadı. Arızalı/asimetrik çalışma için bu indirgeme yeniden değerlendirilmelidir.
- Güç elektroniği PWM/ZOH ve elektromanyetik geçişler 50 Hz statik port hesabına dahil değil. Bunları temsil ettiği iddia edilmiyor.

## 5. Per-unit baz denetimi

\[
I_b=\frac{S_b}{\sqrt3 V_{LL,b}},\quad Z_b=\frac{V_{LL,b}^2}{S_b}.
\]

| Sistem | Sbase | VLL base | Ibase RMS | Zbase |
|---|---:|---:|---:|---:|
| Converter / seçilen ortak LV portu | 2 MVA | 690 V | 1673.479041 A | 0.23805 ohm |
| Fiziksel grid ölçümleri | 2 MVA | 20 kV | 57.735027 A | 200 ohm |
| Trafoyla uyumlu HV referans bazı | 2 MVA | 20.4 kV | 56.602968 A | 208.08 ohm |

RMS fazörlerin pu büyüklüğü, kontrolörün kullanılan tepe-dq bazının pu büyüklüğüyle aynıdır. Ölçüm ölçekleri `sqrt(3/2)/VLL` ve `sqrt(3/2)*VLL/Sn`: faz gerilim/akım tepe değerlerini bazlarına böler. Bu yüzden limiterin 1.2 pu değeri fiziksel olarak **yaklaşık 2008.175 A faz RMS** akımına karşılık gelir; tepe akım bunun sqrt(2) katıdır. 1.2 pu, 100 MVA sistem bazına ait değildir.

Genel dönüşüm:

\[
Z_{pu,new}=Z_{pu,old}\frac{S_{new}}{S_{old}}
\left(\frac{V_{old}}{V_{new}}\right)^2.
\]

Örneğin önceki hat+kaynak `0.00700007438+j0.02706043346` pu, 20 kV bazındaydı. 20.4 kV bazına geçerken `(20/20.4)^2` ile çarpılır. Trafo sargı parametreleri zaten kendi nominal winding bazındadır; aynı ölçek iki kez uygulanmaz.

Grid ölçümlerinde Vpu_20kV=1.02 Vpu_20.4kV, Ipu_20kV=Ipu_20.4kV/1.02; çarpımları aynı kalır. P/Q pu, her iki tarafta da 2 MVA bazındadır. `powergui/Pbase=100e6` veya VI Measurement bloklarının pu çıkışı kapalıyken saklanan Pbase alanları, fiziksel SI ağı veya açık kontrolör ölçeklemesini 100 MVA bazına dönüştürmez. Kaynak kısa-devre gücü 100 MVA da bir kontrolör güç bazı değildir.

## 6. Şöntleri koruyan Thevenin indirgemesi

I_conv = I_trafo_primary + Yc Vc. Limiter, trafo akımını değil **converter-side akımı** izleyen akım döngüsünün referansını sınırlar. Bu yüzden Yc'yi tamamen silmek yanlış akımı verir. Vc'nin port akımına paralel Yc dahil edilerek tam bir Thevenin port eşdeğeri oluşturulabilir; bunun bedeli yalnız Z'nin değil kaynak geriliminin de dönüşmesidir.

Aşağıdaki ardışık indirgeme, 4 iç düğümlü admittans çözümüyle bağımsız olarak karşılaştırıldı:

```text
Za = 1 / (1/Zsource + Yline/2 + Yfault)
Ua = Ubus / (1 + Zsource*(Yline/2 + Yfault))
Zb = 1 / (Yline/2 + 1/(Zline+Za))
Ub = Ua / (1 + (Yline/2)*(Zline+Za))
Zm = 1 / (Ym + 1/(Zt2+Zbreaker+Zb))
Um = Ub / (1 + Ym*(Zt2+Zbreaker+Zb))
Zp = Zt1 + Zm
Zeq = Zp / (1 + Yc*Zp)
Uth = Um / (1 + Yc*Zp)
```

Burada Ubus=0.980392156863∠0 pu, fiziksel HV kaynak LV'ye referanslandıktan sonra seçilen açısal sıfırdır.

\[
\underline U_{th}=0.983679763315-j0.000921202906
=0.983680194662\angle(-0.053656713408^\circ).
\]

\[
\underline I_{conv}=\frac{\underline V_c-\underline U_{th}}
{0.016884633807+j0.066199725744}.
\]

Bu port bağıntısı, tutulan lineer pozitif-sıra ağ için tamdır; **Uth yerine 0.980392 veya 1 yazılamaz**. Sadece seri elemanları toplarsak `0.016776311399+j0.066009643850` elde edilir; fakat bu toplama aynı Uth'yi üretmez. Nominalde bu naif yaklaşımın akım büyüklüğü tesadüfen yakın olsa da **kompleks akım hatası yaklaşık 5.02%**; kayıtlı 0.1 pu yükte yaklaşık 41.64%'tür. Şöntler özellikle reaktif akım ve düşük yükte önemlidir.

## 7. Akım, güç ve limiter sınırı denklemleri

### 7.1 Akım ve açı tanımı

`gamma = arg(Vc)-arg(Uth)` ve `U=abs(Uth)` seçilince:

\[
I_{req}(\gamma)=\frac{E e^{j\gamma}-U}{R+jX},\qquad
|I_{req}|=\frac{\sqrt{E^2+U^2-2EU\cos\gamma}}{\sqrt{R^2+X^2}}.
\]

Bu dönüş yalnız referans eksenini döndürür; fiziksel akım büyüklüğü değişmez. E sabit ve ölçüm kazancı bir kabul edilirse:

\[
\gamma_{lim}=\arccos\left(\frac{E^2+U^2-I_{max}^2(R^2+X^2)}{2EU}\right).
\]

Nominal hesapta arccos argümanı **0.996699705113**, dolayısıyla gerçek kök vardır; ilk pozitif kök **4.65622027424°**. Sabit E için negatif simetrik kök ve periyodik kopyalar da vardır; burada ilk pozitif dal inceleniyor.

Genel fizibilite: Imax < abs(E−U)/abs(Z) ise hiç inaktif açı yoktur; Imax > (E+U)/abs(Z) ise tüm açılar limitin altındadır; eşitlikler teğet sınır durumudur. Sınır dışı arccos argümanı sessizce [-1,1]'e kırpılmamalıdır.

Q–V korunursa E=E(gamma) olduğundan tek bir nominal E ile arccos sonucu tüm ailenin sınırı değildir. `E=1−0.1Qg` ve `abs(Iconv)=1.2` birlikte çözülür:

\[
E_{lim}=1.000990559072,\quad
\gamma_{lim}=4.628267728730^\circ,\quad
P_{g,lim}=1.184595576796\ \mathrm{pu}.
\]

Fiziksel bara referanslı aynı eşik 4.57461101532°'dir. Bu bir **akım limiteri aktivasyon eşiği**; kararlılık veya contraction sınırı değildir.

### 7.2 Güç ölçüm yeri

`measurement_transformations/Power_Calculation_(Positive_Sequence)` şu sinyalleri kullanıyor:

```text
p_grid = v_grid_dq+_d*i_grid_dq+_d + v_grid_dq+_q*i_grid_dq+_q
q_grid = v_grid_dq+_q*i_grid_dq+_d - v_grid_dq+_d*i_grid_dq+_q
```

Gerilim grid terminalinde/kesici sonrasında; akım trafo çıkışında/kesici öncesinde ölçülüyor. Normal kapalı kesicide arada başka dal olmadığından aynı akımdır. VSM swing denklemi **bu p_grid** sinyalini kullanır; Vc-port gücünü veya köprü gücünü kullanmaz.

Port gücünün doğru genel formülü:

\[
P_{port}=\Re\{E e^{j\gamma}I_{conv}^*\}
=\frac{R(E^2-EU\cos\gamma)+XEU\sin\gamma}{R^2+X^2}.
\]

Fakat swing Pe için bunu doğrudan kullanmak uygun değildir. Ağdan grid düğümü ve akımı geri oluşturulmalıdır. Fiziksel bara referansında `v=E exp(j gamma_bus)`, `u=0.980392156863` için bu çalışmada bulunan katsayılar:

```text
vg = (0.394534462826 - j0.001701124220)*v
   + (0.605715891198 + j0.001635965002)*u
ig = (3.615958274343 - j14.227431855063)*v
   + (-3.616905793832 + j14.234702848129)*u
Sg = vg*conj(ig)
Pe = real(Sg),   Qg = imag(Sg)   [birim ölçüm yaklaşımı]
```

Katsayıların tam sayısal halleri MAT/JSON özetinde bulunuyor. Bunlar yeni dinamik durum eklemez. E'nin cebirsel Q–V bağıntısı bu Qg ile çözülür.

Nominalde Pport=1.012440266632 pu iken Pg=1.000000000002 pu. Aradaki yaklaşık **0.01244 pu**; kapasitör damping kaybı, trafo/core kaybı ve grid ölçümüne kadar olan kesici kaybıdır. Hat/kaynak kayıpları ölçümün sonrasında kalır; Pg'den ayrıca çıkartılarak swing'e sonsuz-baraya teslim edilen güç konmamalıdır.

## 8. Ölçüm filtrelerinin statik etkisi: ek hassasiyet kontrolü

Modelde ölçüm LPF'si 2500 Hz, alpha-beta DC giderme HPF'si 2 Hz, Ts=1e-4 s ve DDSRF pozitif sıra çıkarımı aktiftir. Dengeli kararlı durumda DDSRF pozitif bileşeninin DC kazancı 1 kabul edilerek, modeldeki ayrık filtre katsayılarından 50 Hz'de:

\[
G_{50}=0.998801953684+j0.019955151101,
\quad |G|=0.999001276645,
\quad \arg G=1.144565088590^\circ.
\]

Burada `alpha=(1−Ts*pi*2500)/(1+Ts*pi*2500)` LPF katsayısı; HPF `butter(1,2*2*Ts,'high')` karşılığıdır. Fiziksel E'ye göre ölçüm `G*Vc` olur. Böylece ideal sıfır PI hatasıyla:

```text
Ephysical = Eref/abs(G)
arg(Vc) = theta_VSM - arg(G)
Pe_meas = abs(G)^2 * real(vg*conj(ig))
Qe_meas = abs(G)^2 * imag(vg*conj(ig))
Iref magnitude = abs(G)*abs(Iconv)
```

Bu, filtre durumlarını modele eklemek değil, yalnız temel frekans statik kazanç/faz düzeltmesidir. Rate-transition zamanlaması, ZOH ve EMT ayrıklaştırma etkilerinin tamamını temsil ettiği iddia edilmiyor.

| Büyüklük | Statik ölçüm düzeltmesiyle |
|---|---:|
| Fiziksel nominal E | 0.999422979702 pu |
| Eref | 0.998424832631 pu |
| Fiziksel grid P / ölçülen P | 1.002000443045 / 1.000000000000 pu |
| Ölçülen Q | 0.015751673690 pu |
| Fiziksel abs(Iconv) / abs(Iref) | 1.015155467085 / 1.014141607611 pu |
| Fiziksel akımda 1.2 referansın karşılığı | 1.201199666161 pu |
| Sabit nominal E ile gamma_lim | 4.655387145755 deg |
| Cebirsel Q–V ile gamma_lim | 4.626919841582 deg |

Son eşiğin **VSG açısı − LV'ye referanslanmış fiziksel bara açısı** karşılığı 5.717828216764°'dir; çünkü fiziksel açıya arg(G) eklenir. Bu nedenle 4.6269°, 5.7178° ve ham HV referanslı theta değerleri aynı koordinat değildir. 30° trafo dönüşü ve filtre fazı sessizce göz ardı edilmemelidir.

## 9. Yapısal ve sayısal tutarlılık kontrolleri

- Düğüm admittans çözümü ile bağımsız ardışık Thevenin indirgemesi uyuştu. Birden fazla E ve açı için en büyük KCL/Thevenin artık büyüklüğü yaklaşık 7.3e-12 pu.
- Nominal yük akışı ve Q–V denklemi birlikte çözüldü; residual yaklaşık 1.8e-12 pu veya altında.
- Nominal abs(Iref)=1.01348 < 1.2; hesaplanan köprü gerilimi 1.00173 pu. Aktif sönümleme eklenmeden önce gerekli ham gerilim referansı yaklaşık **1.20413 pu**, yani 1.3'ün altında.
- Akım aktivasyon eşiğinde gerekli ham gerilim referansı yaklaşık **1.24385 pu**, statik ölçüm düzeltmesiyle 1.24474 pu; dolayısıyla bu ideal normal-takip hesabında gerilim limiteri akım limiterinden önce devreye girmiyor. Bu kontrol, limiter-inaktif taraftaki varsayımı destekler; genel EMT geçiş garantisi değildir.
- E=1 varsayımı sonuç olarak nominale yakındır, fakat gerekçesi Qstar=0 değil; çözülmüş Qg ve droop ilişkisidir.
- Kayıtlı P_set seviyelerinde aynı nominal-fazör yaklaşımının verdiği değerler:

| Pstar pu | E pu | Qg pu | gamma_bus deg | abs(Iconv) pu |
|---:|---:|---:|---:|---:|
| 0.1 | 0.989577274 | 0.104227261 | 0.283066330 | 0.121052355 |
| 0.35 | 0.992151274 | 0.078487256 | 1.279611955 | 0.358598699 |
| 0.4 | 0.992670198 | 0.073298024 | 1.478365102 | 0.408511771 |
| 1.0, analitik nominal | 0.999003208 | 0.009967915 | 3.849790563 | 1.013477853 |

Bu tablo EMT ölçümü değil, model tabanlı öngörüdür. Kayıtlı modelde Pstar=1 noktasının simüle edilmiş olduğu söylenemez.

## 10. Eski 20.7395° sonucunun yorumu

Eski değer E=U=1, X=0.3 ve Imax=1.2 tanımında matematiksel olarak doğrudur. Ancak bu X fiziksel ağın doğrulanmış eşdeğeri değildi; ayrıca aynı akım, güç ölçüm yeri ve trafo/şönt etkileri korunmamıştı.

Yeni eşdeğerde abs(Z)=0.068319 ve limit eşiği kullanılan E yaklaşımına göre yaklaşık 4.63–4.66° (Thevenin referanslı). Dolayısıyla **eski 20.74° fiziksel modele yakın değildir**. Eşlenmiş fiziksel/Thevenin açı koordinatında limiterin başlamasını yaklaşık 4.5 kat daha ileri öngörür; limiter-inaktif açı aralığı açısından iyimser kalır. Bu, eski contraction aralığının fiziksel güvenlik marjı olduğu anlamına gelmez; iki ayrı elektriksel varsayımın sonucudur. Ham VSM/HV açısını doğrudan 4.63° ile kıyaslamak da doğru değildir.

## 11. İki durum korunarak önerilen minimum revizyon

Karar: **elektriksel cebirsel modeli yapısal olarak düzeltmek, sonra analizi yeniden çalıştırmak**. Sadece Xeq değerini değiştirip eski sonuçları fiziksel sonuç olarak sunmak uygun değil.

1. Dinamik durumlar `[gamma,nu]` olarak kalabilir. Seçilen açı referansı açıkça tanımlanmalı; fiziksel Vc açısı ile VSG theta arasındaki sabit trafo/ölçüm ofsetleri kaydedilmeli.
2. Normal çalışmada Vc genliği `E=1−0.1Qg` cebirsel koşulundan, grid P/Q ise yukarıdaki lineer ağ katsayılarından hesaplanmalı. İstenirse önce açıkça etiketlenmiş sabit nominal E yaklaşımı kullanılabilir; MAT'taki üst seviye E/U/Z/gamma_lim değerleri bu yaklaşım için birbirleriyle uyumludur.
3. Limiter talebi converter akımı üzerinden bulunmalı: `Ireq=(Vdes−Uth)/Zeq`. Grid/trafo akımıyla karıştırılmamalı.
4. Aktif limiteri daha sonra aynı iki durumlu ideal modelde korumak istenirse, pasif ağın Kirchhoff bağıntılarını koruyan minimal aday: `Ic_lim=sat(Ireq)`, `Vc_actual=Uth+Zeq*Ic_lim`, ardından grid düğümleri ve Pe/Qe'yi **Vc_actual** üzerinden geri kurmak. Sınırlanan akıma rağmen Vc_actual=Vdes dayatmak iki uçlu pasif ağda genel olarak tutarsızdır. Bu aday yeni bir ideal cebirsel kontrol kapanışıdır; mevcut görevde uygulanmadı/doğrulanmadı.
5. Limiter aktif olduğunda fiziksel Vc'nin açısı artık yalnız VSG açısına eşit olmayabilir. E ve gamma'yı o dalda hedef gerilim mi yoksa gerçekleşen gerilim mi olarak kullandığımız açıkça seçilmelidir. Gerçek PI/freeze EMT modelinin birebir eşdeğerliği iddia edilmemeli.
6. Swing Pe, `Re(Vc*Ic*)` yerine modeldeki **grid tarafı güç ölçümüne** karşılık gelmeli. Özellikle limiterde `Pe=rho*eski_sinus` ilişkisi otomatik olarak korunamaz.

İlk sonraki çalışma, bu cebirsel kapanışı ve açı tanımını sabitleyip yeni denge ailesini üretmek olmalıdır. Bu raporda yeni contraction bölgeleri hesaplanmadı ve önceki bölgeler fiziksel model için yeniden etiketlenmedi.

## 12. Çıktılar ve tekrar üretim

- `gfm_E_Zeq_identification_summary.mat`: istenen E_pu, U1_pu, Req_pu, Xeq_pu, Zeq_mag_pu, Zeq_angle_deg, Imax_pu, gamma_lim_deg ve ortak baz alanlarını üst seviyede içerir. **U1_pu burada Thevenin kaynak büyüklüğüdür**; fiziksel bara ayrıca `U1_physical_bus_referred_pu` alanındadır. Üst seviye gamma_lim sabit nominal E'ye aittir; değişken E/Q–V sonucu ayrı `threshold_algebraic_QV` yapısındadır.
- `gfm_E_Zeq_identification_summary.json`: aynı sayısal içerik, okunabilir biçimde.
- `gfm_E_Zeq_identification_console.txt`: eski/yeni karşılaştırması ve karar.
- `diagnostics/gfm_E_Zeq_identification.py`: parametreler, düğüm çözümü, bağımsız Thevenin denetimi ve sayısal kontroller. Mevcut MATLAB analizlerine yazmaz.

Korunan dosyaların SHA-256 değerleri:

```text
baseline_gfm_contraction.m
c532dca6af2188fa0bd8fd8fffca74e0d02d6f345df0ebd368fb9a0401b128d0
gfm_current_limiter_analysis.m
8fecab20a0ead40f07d10bcba4315ab85418c01a662f246eaa2fb7e404e312ed
```
