# GFMVSG1: gerçek akım sınırlayıcısının indirgenmesi için teknik inceleme

İnceleme tarihi: 30 Eylül 2026. Bu çalışma dosyalardaki blok, parametre ve bağlantı incelemesidir; yeni bir Simulink çalışma noktası çözümü veya simülasyon doğrulaması değildir. `baseline_gfm_contraction.m` ve `gfm_current_limiter_analysis.m` değiştirilmedi; yeni contraction bölgesi hesaplanmadı.

**Ana sonuç:** gerçek yapı, akım referansının radyal doyurulması ve PI integratörlerinin eksen bazında koşullu dondurulmasıdır. İki durumlu ideal akım kırpma modeli bunun eşdeğeri değildir. Altı durumlu C adayı kontrolör hafızasını korumak için uygun bir başlangıçtır; fakat gerçek Simulink sisteminin tam, kapalı ve gecikmesiz modeli değildir. Elektrik devresi, ölçüm filtreleri ve dış kontrolün nasıl indirgeneceği ayrıca belirlenmelidir.

## İncelenen kaynaklar ve doğrulama kapsamı

Ana model: `model/GFMVSG1.slx`.

Kütüphane kökü: `model/tdmsquared-master/tdmsquared-master/Library`.

- `Distributed Generation Units/Converter/Grid_Forming_Converter/Subsys_GFormConv_InnerCtrl_cascadedPI.slx`: `v_cap_Controller` (SID 91), `i_conv_Controller` (SID 21), freeze sistemleri 6708/6735, force-freeze sistemi 65, akım sınırlayıcı sistemi 661.
- Aynı klasörde `Subsys_GFormConv_Control.slx`: çalışma modu, sanal endüktans ve aktif sönümleme bağlantıları.
- Aynı klasörde `Subsys_GFormConv_OuterCtrl_VSM.slx`: frekans, açı ve Q–V dış kontrolü.
- `Shared Functions/PI_freeze.slx`, `Shared Functions/dq_Limiter.slx`.
- `Distributed Generation Units/Converter/Converter_Shared_Subsystems/Subsys_Conv_Electrical_Part_SPS.slx` ve `Subsys_Conv_Ctrl_measurement_transformations.slx`.

İkinci `code - Simulation/.../Version 2_3/.../Library` ağacındaki iç kontrol ve PI_freeze sistem XML'leri karşılaştırıldı: bulunan farklar pencerenin açık olması ve çizim sırası gibi görünüm bilgileridir; incelenen denklemleri değiştirmiyor. dq_Limiter ve ana kontrol dosyaları iki ağaçta birebir aynıdır. Bu rapor, MATLAB çalışma anındaki path seçiminin bütünü için bir doğrulama iddiası taşımaz.

## A. Modelden çıkarılan kontrolör denklemleri

### A1. Birimler, eksenler ve gösterim

Gerilim, akım, güç ve frekans kontrol sinyalleri pu'dur. Frekans değişkeni `nu = f_pu`'dur; iç kontrolün çapraz bağ terimleri Hz ile çarpılmıyor. `f_n` kazancı ayrı izleme dalındadır. Nominal frekans 50 Hz, güç bazı 2 MVA, dönüştürücü hatlar arası gerilim bazı 690 V, şebeke bazı 20 kV'dur. Ölçüm dönüşümünde gerilim ölçeği `sqrt(3/2)/V_LL`, akım ölçeği `sqrt(3/2)*V_LL/S_n` kullanılır.

PI girişleri çalışma modu tarafından seçilen theta açısındaki **pozitif sıralı dq+ ölçümleridir**. Normal GFM modunda açı VSM'den gelir. Aşağıdaki `v_c`, `i_c`, `i_f` sırasıyla `v_cap_dq+`, `i_conv_dq+`, `i_filter_grid_dq+` ölçümleridir. `i_f` ile trafo sonrası `i_grid_dq+` aynı sinyal değildir. Aktif sönümleme ayrı olarak `i_conv_dq0` kullanır; geçişlerde bu sinyali dq+ ölçümüyle özdeş saymak bir yaklaşımdır.

`xi_vd, xi_vq, xi_id, xi_iq`, PI integratörlerinin çıkış katkıları olarak tanımlanmıştır; dolayısıyla türev/güncellemede Ki bulunur. Başlangıç/reset değeri sıfırdır; bu, yüklü çalışma noktasındaki denge değerlerinin sıfır olduğu anlamına gelmez.

### A2. Gerilim PI → ham akım referansı

`e_vd = v_d* - v_cd`, `e_vq = v_q* - v_cq`.

```text
r_d = 0.8 e_vd + xi_vd + i_fd - nu (0.04569) v_cq
r_q = 0.8 e_vq + xi_vq + i_fq + nu (0.04569) v_cd
```

Burada `r = i_ref_raw`. Akım feedforward etkin (`K_ff_I_grid=1`), referans türevi command-feedforward kapalıdır. `0.04569 = C_filter_pu(2)` kapasitif çapraz bağ katsayısıdır. Ölçülen gerilim doğrudan PI hatasına; filtre çıkış akımı feedforward dalına girer.

`OuterCtrl_Lvirt=0` olduğundan normal GFM modunda `v_d*=V_outer`, `v_q*=0`. Ancak gerçek modelde `V_outer` sabit E=1 olarak alınamaz: `V_outer = v_set + y_Q`; `y_Q`, `0.1(q_set-q_grid)` girişinin T_Q=0.2 s dış kontrol filtresinden geçen değeridir. Gerçek ayrık filtre transferi `(1-alpha_LPF_Q)(1+z^-1)/(2(1-alpha_LPF_Q z^-1))` biçimindedir. Sürekli zaman `0.2 dy_Q/dt = -y_Q + 0.1(q_set-Qe)` ifadesi bu bloğun indirgeme karşılığıdır, ayrık bloğun birebir denklemi değildir.

### A3. Akım PI → ham gerilim referansı

`s = i_ref_lim`, `e_id=s_d-i_cd`, `e_iq=s_q-i_cq`.

```text
u_d = 3.3 e_id + xi_id + v_cd - nu (0.05097) i_cq
u_q = 3.3 e_iq + xi_iq + v_cq + nu (0.05097) i_cd
```

Gerilim feedforward etkin (`K_ff_V_cap=1`); akım referansı türevi command-feedforward kapalıdır. `0.05097 = L_filter_conv_pu(2)`. Akım ölçümü PI hatasında ve çapraz bağlarda kullanılır. İncelenen toplama bloklarında ayrıca bir `R_filter*i_c` feedforward terimi yoktur; fiziksel direnç elektrik devresinde bulunur.

Gerilim sınırlayıcısından SONRA aktif sönümleme uygulanır:

```text
w = sat_radial(u, 1.3)
v_cmd,dq = w - 0.2 i_conv,dq0
```

Ardından dq0→abc dönüşümü ve elektriksel dönüştürücü yolu gelir. **1.3 sınırı `u→w` noktasındadır; `v_cmd` büyüklüğünün 1.3'ü aşmayacağını garanti etmez.** Aktif sönümleme sonrası çıkış bu gerilim sınırlayıcısının freeze algılayıcısına geri beslenmiyor. `0.2` aktif sönümleme direncidir; Unified-paper SatLim anti-windup `Z_vi` parametresi değildir.

### A4. Dış VSM ile baseline ayrımı

Normal etkin çalışma için blokların sürekli zaman karşılığı:

```text
dnu/dt = [(p_set-p_grid)/max(nu,0.5)
          + K(f_set_pu-nu) - Kd(nu-f_PLL)]/(2H)
dtheta/dt = 2*pi*50*nu
```

VSM frekans integratörü gerçekte ayrık olup mod kontrollü reset ve harici PLL başlangıç değeri içerir. Kaydedilmiş Simulink örneğinde H=0.1, K=20, Kd=40; doğrulanmış analitik baseline'da H=3.4, K=20, Kd=5 kullanılıyor. Bunlar farklı tuning'lerdir.

`f_PLL=f_set_pu=1` varsayımı ve `nu>0.5` altında sönümleme `(K+Kd)(nu-1)` olur. Gerçek modelde PLL'nin varlığı yalnız başına GFL kontrolü anlamına gelmez; burada VSM sönümleme ve senkronizasyon yollarında da kullanılır. PLL'yi sabit şebeke frekansına indirgemek açık bir yaklaşım olarak belirtilmelidir.

## B. Tam limiter ve freeze mantığı

### B1. Radyal doyum

Her iki limiter için:

```text
sat_radial(a,c) = a                         if ||a|| <= c
                = c*a/||a||                if ||a|| > c
rho(a,c)        = 1                         if ||a|| = 0
                = min(1,c/||a||)           otherwise
```

Akım için `c=1.2`, gerilim için `c=1.3`. Seçili dq limiter `Constant ratio d/q` seçeneğidir. Büyük norm halinde açı `atan2` ile bulunup sin/cos ile sınır çemberine taşınır; yukarıdaki formülle eşdeğerdir. d/q önceliklendirmesi seçili değildir. Sınırlanmış vektörün kendisine bu blokta bir örnek gecikmesi konmamıştır; aktivasyon/freeze yolunda gecikme vardır. Sınırda aktivasyon testi `>` olduğundan eşitlikte limiter aktif bayrağı false'dur.

### B2. Ayrık freeze denklemleri

`a_I[k]=(||r[k]||>1.2)`, `a_V[k]=(||u[k]||>1.3)` tanımlansın. d veya q ekseni için:

```text
F_va[k] = F_force[k] OR
          (a_I[k-1] AND (sign(r_a[k-1]) == sign(e_va[k-1])))
F_ia[k] = F_force[k] OR
          (a_V[k-1] AND (sign(u_a[k-1]) == sign(e_ia[k-1])))
```

Karşılaştırma ve aktivasyon paralel yollarda birer örnek gecikir; bunlar iki örneklik seri gecikme değildir. `F_force=1` mod -1 ve 1 için; `F_force=0` mod 0 ve 2 için geçerlidir. Bu eşleme switch bağlantılarından doğrulandı.

Reset olayı yokken Ts=1e-4 s ve Forward Euler ile:

```text
xi_va[k+1] = xi_va[k] + Ts*0.5 *(1-F_va[k])*e_va[k]
xi_ia[k+1] = xi_ia[k] + Ts*16.5*(1-F_ia[k])*e_ia[k]
```

Force-freeze reset sinyalinin yükselen kenarında ilgili PI integratörleri sıfıra resetlenir. Freeze sırasında proporcional katkı `Kp*e` çalışmaya devam eder: PI_freeze harici bypass yolu bunu korur. İç PID bloğunda `AntiWindupMode=none`, tracking ve dahili çıkış/integratör sınırı kapalıdır.

Normal modda limiter inaktifse integrasyon sürer; aktif ve işaretler farklıysa integrasyon sürer; aktif ve işaretler aynıysa ilgili eksen donar. `sign(0)==sign(0)` doğrudur; buna rağmen hata sıfırsa artım zaten sıfırdır. Sonuç eksen bazındadır; radyal normun tek başına türev testiyle değiştirilmemelidir.

Gecikme ve örnekleme ihmal edilirse `dot(xi_va)=0.5(1-F_va)e_va`, `dot(xi_ia)=16.5(1-F_ia)e_ia` yazılabilir. Bu, gecikmesiz hibrit sürekli-zaman yaklaşımıdır; tam Simulink uygulaması değildir.

**Gerilim limiterinin ilgili çalışma noktalarında inaktif olduğu doğrulanmış değildir.** Sadece Imax ve gamma bilgisi bunu belirlemeye yetmez; `u_d,u_q` ve integral durumları gereklidir. Dolayısıyla current-PI durumlarını sırf 1.3 sınırını büyük kabul ederek atamayız.

## C. Önerilen minimum durum vektörü

| Aday | Hafıza ve Markov özelliği | Freeze ve Jacobian açısından karar |
|---|---|---|
| A: `[gamma,nu]` | Gerilim ve akım PI hafızasının tamamını atar. Gerçek kontrolör için kapalı değildir. | İdeal statik sınırlayıcı modeli olabilir; gerçek freeze'i temsil etmez. |
| B: `[gamma,nu,xi_vd,xi_vq]` | Gerilim PI hafızasını korur, akım PI hafızasını atar. | Yalnızca ayrıca doğrulanmış ideal/hızlı akım döngüsü indirgemesiyle kapatılabilir. Gerilim limiteri/current-PI freeze'i genel olarak temsil edilemez. |
| C: `[gamma,nu,xi_vd,xi_vq,xi_id,xi_iq]` | İki PI çiftinin hafızasını korur; fakat elektrik devresi, dış Q filtresi, PLL/ölçüm filtreleri ve örnek gecikmeleri hâlâ eksiktir. | Açıkça tanımlanmış cebirsel ağ ve dış kontrol varsayımlarıyla DAE çekirdeği olarak önerilir. Tek başına tam Simulink modeli diye sunulmamalıdır. |

**Öneri:** C'yi altı durumlu kontrolör çekirdeği olarak kullanmak; önce elektrik dinamiklerini koruyan bir karşılaştırma modeliyle bunun cebirsel indirgemesini sınamak. Gerçek dış Q–V kontrolü korunursa en az `y_Q` da eklenir. Elektrik dinamikleri tutulursa dönüştürücü RL akımı, şönt RC kolunun kapasitör gerilimi ve trafo/şebeke akımları gibi bağımsız enerji depolama durumları gerekir. Tam ağın hat kapasiteleri, trafo mıknatıslanması ve DC tarafı da seçilen doğruluk seviyesine göre ek durum getirir.

Tam ayrık kontrol için aktivasyon ve dört işaret karşılaştırmasının önceki örnek değerleri, reset geçmişi ve ölçüm/filtre durumları da durum belleğinin parçasıdır. Bu nedenle eldeki inceleme, tam model için kesin bir “minimum altı durum” sonucu vermiyor. Yaklaşım kapsamı seçilmeden tek bir minimum boyut ileri sürmek doğru değildir.

## D. Cebirsel değişkenler ve elektriksel kapanış

Quasi-steady dengeli ağ yaklaşımı kabul edilirse bağımsız cebirsel bilinmeyenler örneğin:

```text
z = [i_conv,dq; i_filter_grid,dq; v_cap_node,dq;
     i_grid,dq; v_grid_terminal,dq; gerekli ara ağ düğümleri]
```

Trafo oranı ve faz kayması, şönt kol ve hat şöntleri tutuluyorsa bunlara ilişkin düğüm/kol bağıntıları da g içine girer. Buradaki düğüm gerilimi ile seri RC kolundaki kapasitörün kendi gerilimi birbirinden ayrılmalıdır.

`e_v,e_i,r,s,u,w,v_cmd,||r||,||u||,rho_I,rho_V,Pe,Qe` hesaplanabilen türetilmiş büyüklüklerdir; çözüm düzeni için z içine konmaları şart değildir. Ölçüm filtreleri korunuyorsa ölçülen dq+ sinyaller anlık fiziksel değişkenlerden cebirsel elde edilemez; filtre durumlarıyla üretilir.

Gerçek güç ölçümü trafo sonrası `v_grid_dq+` ve `i_grid_dq+` kullanır:

```text
Pe = v_gd*i_gd + v_gq*i_gq
Qe = v_gq*i_gd - v_gd*i_gq
```

Fiziksel yol: sınırlanmış gerilim referansı → aktif sönümleme → ters Park/dönüştürücü → dönüştürücü tarafı RL → şönt seri-RC filtre düğümü → seçili trafo → kesici/ölçüm → PI-section hat → sonlu kaynak empedansı → ideal kaynak. Akım referansının 1.2 pu ile sınırlandırılması, geçiş anlarında gerçek ölçülen dönüştürücü akımının aynı değerde sert kesildiği anlamına gelmez.

Elektriksel parametreler:

- Dönüştürücü RL: R=0.001 pu, X=0.05097 pu.
- Şönt seri RC: R=0.06875 pu, B_C=0.04569 pu. Kapasitans bloğunda `C=B_C*S_n/(2*pi*50*V_LL^2)` kullanılıyor.
- Trafo iki sargı: her biri R=0.005 pu, X=0.02 pu; Rm=Lm=500 pu. Bağlantı Yg/Delta(D1), nominal gerilimler 690 V / **1.02*20 kV**. Faz kayması ve oran cebirsel indirgemede hesaba katılmalıdır.
- `L_filter_grid_pu=[0.01,0.04]` bloğu **No_Transformer** alternatifindedir. Seçili trafoya ek bir seri bobin olarak toplanmamalıdır.
- Hat: 2 km, pozitif sıra R=0.501 ohm/km, L=2.279 mH/km, C=151.2 nF/km. Hat şöntleri saf seri toplamda kaybolur.
- Kaynak aktif kısa-devre parametreleri: Ssc=100 MVA, X/R=10; 20 kV'da R=0.398014876 ohm, L=12.669206991 mH. Alternatif doğrudan R/L alanları aktif kaynak ayarı değildir.
- Yalnızca hat+kaynak seri eşdeğeri: Rg=1.400014876 ohm, Lg=17.227206991 mH; 2 MVA/20 kV bazında Rg=0.00700007438 pu, Xg=0.02706043346 pu.

Bu bilgiler fiziksel ağ indirgemesini kurmak içindir; `Xeq=0.3` yerine doğrudan kullanılacak tek bir eşdeğer belirlemez.

## E. Contraction için gereken yaklaşımlar

Önerilen sürekli-zaman yaklaşımının doğal biçimi:

```text
xdot = f(x,z),       0 = g(x,z)
Jred = fx - fz*(gz\gx)
```

Son ifade yalnızca ilgili çalışma dalında g_z terslenebilir ve cebirsel çözüm yerel olarak tek olduğunda geçerlidir. Sayısal uygulamada açık matris tersi yerine lineer çözüm kullanılmalıdır. Elektrik/ölçüm dinamikleri tutuluyorsa bunları z'ye zorla taşımamak gerekir.

Radyal doyumun türevi inaktif bölgede I, aktif bölgede `c/||a||*(I-aa'/||a||^2)` olur; sınırda iki taraftaki türevler aynı değildir. Sign/freeze geçişleri de düzgün değildir. Dolayısıyla tek bir düzgün Jacobian'ın tüm çalışma alanını açıklaması beklenmemelidir.

İlk tercih parça parça mod analizi ve geçişlerin ayrıca kontrolüdür. Mod içi negatif simetrik Jacobian tek başına tüm hibrit sistemin contraction kanıtı değildir; ortak metrik, uygun geçiş koşulları ve gerekiyorsa saltation/reset etkileri gerekir. Exact ayrık sistemde önce gecikme belleğiyle kapalı güncelleme haritası kurulur; diferansiyel harita için `A_k' M_next A_k - M_current` türü ayrık contraction koşulu incelenir. Sürekli-zaman `Sym(F)` testi doğrudan bu haritanın yerine geçmez.

Yumuşatma ancak parametresi açıklanan ve orijinal anahtarlamayla karşılaştırılan bir sayısal yaklaşım olarak önerilebilir; modelden çıkarılmış gerçek limiter diye sunulmamalıdır.

Ek zorunlu seçimler: dengeli pozitif sıra varsayımı, ölçüm filtrelerinin korunması/indirgenmesi, PLL frekansının sabitlenmesi, Q–V dinamiğinin korunması, DC/dönüştürücü gerilim üretiminin idealleştirilmesi ve elektrik devresinin quasi-steady kabulünün doğrulanması.

**Zaman birimi:** Ts saniyedir. Fiziksel şebeke için `dot(gamma)=2*pi*50*(nu-1)`; mevcut baseline ise açıkça belirtilmiş normalize analitik zamanda `dot(gamma)=nu-1` kullanıyor. Yalnızca açı denklemindeki katsayıyı değiştirmek tutarlı zaman dönüşümü değildir. PI kazançları, dış dinamikler ve Ts ortak zaman ölçeğine dönüştürülmeden bu modeller birleştirilmemelidir.

## F. Xeq=0.3 limiter-free baseline tutarlılığı

Baseline dosyasında E=U1=1, Xeq=0.3 doğrudan atanıyor. İncelenen satırlarda bunun filtre, trafo ve hat verilerinden çıkarıldığı bir fiziksel eşdeğerleme yok. Dolayısıyla 0.3'ün bu devrede hangi bileşenlerin toplamı olduğunu doğrulayamadık; şu aşamada bağımsız indirgenmiş model parametresidir.

Baseline ilişkisi, iki ideal gerilim kaynağı arasındaki kayıpsız saf reaktans için elde edilir:

```text
i = (E*exp(j*gamma)-U1)/(j*Xeq)
Pe = real(E*exp(j*gamma)*conj(i))
   = E*U1/Xeq*sin(gamma)
```

Dirençli tek seri dalda bile gönderici güç bağıntısı
`[R*(E^2-E*U1*cos(gamma))+X*E*U1*sin(gamma)]/(R^2+X^2)` olur. Gerçek ağda buna trafo oranı/fazı, şönt dallar, ölçüm noktası ve dış gerilim kontrolü de eklenir. Kaynak tarafı güç ile p_grid ölçümü dirençli ağda aynı değildir.

Bu nedenle fiziksel ağı koruyan quasi-steady modelin yalnız limiteri kapatarak baseline sinüs bağıntısını aynen geri verdiği **henüz gösterilmiş değildir**. Bu kontrol şu anda geçilmiş sayılmaz. Baseline'a tam dönüş istenirse ayrıca R=0, uygun referans açı tanımı, sabit E, ideal gerilim takibi, ihmal edilmiş şöntler ve eşdeğer X=0.3 gibi açık varsayımlar gerekir. Bu durumda oluşan model baseline ile uyumlu bir idealizasyon olur; gerçek devrenin birebir indirgemesi olarak adlandırılmamalıdır.

## G. Dört noktalı sayısal kontrol

Gerçek PI kontrolöründe gamma tek başına çalışma noktasını belirlemez. `xi_v,xi_i`, v_cap, iki ayrı akım ölçümü, Q/V setpoint ve dış filtre durumu, PLL/ölçüm durumları ve gecikmiş freeze bilgileri gerekir. İncelenen dört açı için tutarlı bir Simulink operating-point kaydı veya bu tam denklemlerle çözülmüş denge mevcut değil. Gerilim limiteri inaktifliği varsayılarak sayı üretmek döngüsel doğrulama olur.

| gamma0 (deg) | ||r|| | ||s|| | rho_I | Freeze_v,d/q | ||u|| | 1.3 aktif? | Freeze_i,d/q | Gerçek Pe | Tam model dengesi? |
|---:|---|---|---|---|---|---|---|---|---|
| 17.457603 | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Doğrulanmadı | Belirlenemedi | Belirlenemedi | Doğrulanmadı |
| 20.5 | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Doğrulanmadı | Belirlenemedi | Belirlenemedi | Doğrulanmadı |
| 25 | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Doğrulanmadı | Belirlenemedi | Belirlenemedi | Doğrulanmadı |
| 65 | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Belirlenemedi | Doğrulanmadı | Belirlenemedi | Belirlenemedi | Doğrulanmadı |

Karşılaştırma amacıyla **yalnız eski ideal E=U1=1, Xeq=0.3 modeli** için hesaplanabilen sayılar aşağıdadır. Bunlar gerçek ham PI akım referansı veya gerçek PI/gerilim limiter doğrulaması değildir.

| gamma0 (deg) | İdeal I_req | min(I_req,1.2) | İdeal rho | Baseline Pstar=sin(gamma0)/0.3 | İdeal rho*Pe |
|---:|---:|---:|---:|---:|---:|
| 17.457603 | 1.011718 | 1.011718 | 1.000000 | 1.000000 | 1.000000 |
| 20.5 | 1.186290 | 1.186290 | 1.000000 | 1.167358 | 1.167358 |
| 25 | 1.442931 | 1.200000 | 0.831641 | 1.408728 | 1.171555 |
| 65 | 3.581997 | 1.200000 | 0.335009 | 3.021026 | 1.012070 |

Bu eski ideal limiter modelinde ilk iki açı eski Pstar değerleriyle dengede kalabilir; son iki açı aynı Pstar ile dengede değildir. Bu sonuç gerçek PI kontrollü Simulink dengeleri hakkında kanıt değildir.

Sayısal tabloyu tamamlamak için her adayda normal mod, tutarlı tuning ve setpoint'lerle denge çözülmeli; r,s,u,w,xi_v,xi_i,ölçümler,Pe ve gecikmiş freeze sinyalleri alınmalıdır. Hatalar, state/update residual'ları ve ağ residual'ları birlikte denetlenmelidir. Sadece açı ve frekansın sabit görünmesi yeterli değildir.

## H. Bir sonraki model için öneri

1. Altı durumlu C kontrolör çekirdeğini, iki radyal sınırlayıcıyı ve eksen bazlı freeze mantığını temel alın. 0.2 aktif sönümlemeyi gerilim doyumundan sonra tutun; anti-windup Zvi eklemeyin.
2. İlk doğrulama modelinde elektrik ve önemli ölçüm/dış kontrol hafızasını koruyun. Kullanılacak H/K/Kd tuning'ini ve fiziksel/normalize zaman ayrımını açıkça tanımlayın.
3. Dört gerçek dengeyi ve limiter durumlarını çözün. Sonra hızlı durumları cebirsel indirgeme yapmanın hem dengeyi hem geçişleri yeterince koruduğunu kontrol edin; g_z düzenliliğini sınayın.
4. Xeq=0.3 baseline'ı fiziksel modelle aynı kabul etmeyin. Baseline-limit testi için gereken ek idealizasyonları açık bir ayrı tutarlılık testi olarak tanımlayın.
5. Ancak bu kapanış ve doğrulamalardan sonra Jred ve hibrit contraction analizine geçin.

Bu incelemede MATLAB analiz dosyalarında değişiklik yapılmadı. Yeni contraction bölgesi, yeni denge iddiası veya doğrulanmamış gerilim-limiter inaktifliği sonucu üretilmedi.
