# 50-Hz GFM/VSM: gerçek zaman türevli makale yöntemi

Aynı `baseline_gfm_contraction.m` dosyası direktiflere göre güncellendi ve
MATLAB R2024b'de çalıştırıldı. Önceki Lyapunov referans dalı ve
omega0*d(beta)/d(gamma) tasarım katsayısı kaldırıldı. Önceki 82 derece
civarındaki contraction aralıkları artık bu sürümün sonucu değildir.
Eski sonuç ve grafik dosyaları güncel sonuçlarla yenilendi.

## Sabit fiziksel parametreler

f0=50 Hz, omega0=314.159265359 rad/s, Mv=0.0042 s²/rad,
mp_pct=1%, mp=3.14159265359 rad/s/pu, D=0.318309886184 s/rad,
D/Mv=75.788068139 1/s. Rg=0.0139, Xg=0.166, E=Vs=1 pu.
Qstar=0 yalnızca metadata; sabit E modeli reaktif güç regülasyonu içermez.
Akım sınırlama yoktur. Her gamma0 için Pstar0=Pe(gamma0) ayrı hesaplanır.

## Türevler ve dönüşüm

x=[gamma; Delta omega], u=Vs. Mevcut dirençli-endüktif güç bağıntısı korunmuştur.
C=E/(Mv*|Zg|), theta=gamma+phig ve w=Delta omega tanımlarıyla:

- A=[0,1; -C*u*sin(theta), -D/Mv].
- B=[0; C*cos(theta)].
- alpha1=A(2,1), alpha2=A(2,2), beta=C*cos(theta).
- beta_dot=-C*sin(theta)*w.
- beta_ddot=-C*cos(theta)*w²-C*sin(theta)*f2(x,u).
- Theta=[-beta,0; -beta_dot,-beta].
- Theta_dot=[-beta_dot,0; -beta_ddot,-beta_dot].

beta, u'ya doğrudan bağlı değildir. Eq.(12) dönüşümü kodda açıkça uygulanır;
GFM'nin A,B yapısıyla beta sıfır değilken T=I çıkar. Atilde ve Btilde
beklenen companion formlarla sayısal olarak karşılaştırılır.

Gerçek bir dengede w=0 ve f2=0 olduğundan beta_dot=beta_ddot=0,
Theta_dot=0, Theta=-beta0*I, M=beta0²*I ve F=A elde edilir.
Bu sıfırlar varsayılmamış, genel türevlerden çıkarılmıştır.
F=(Theta_dot+Theta*Atilde)/Theta ve Sym(F)=(F+F')/2 doğrudan hesaplanır.
Eq.(15c) biçimindeki analitik matris, yarım-toplam konvansiyonuyla karşılaştırılır.

## Neden Case I yok?

k=Pe'(gamma0)/Mv ve b=D/Mv için

Sym(F) = [0, (1-k)/2; (1-k)/2, -b].

Dolayısıyla det(Sym(F))=-(1-k)²/4 <= 0. Bu metrik kuruluşu ile hiçbir
çalışma noktasında sıkı negatif tanımlılık sağlanamaz. Contraction aralıkları
hem radyan hem derece cinsinden BOŞTUR. Sertifika dalı bu teşhisle durdurulur;
başka bir metrikle değiştirilmez. Geçerli kalan matris teşhisleri ve mevcut
case kategorilerinin yörünge simülasyonları devam eder.

Bu, başka bir metrikle contraction bulunamayacağı veya tüm sistemin kararsız
olduğu anlamına gelmez. Sadece direktifteki dönüşümün mevcut dönen referans
çerçevesindeki iki durumlu GFM modeline bu uygulaması sertifika vermez.

beta=0 noktaları yaklaşık 4.786492591 ve -175.213507409 derecedir.
Buralarda dönüşüm tekildir; epsilon eklenmeden ilgili değerlendirme atlanır.
Yuvarlama seviyesindeki analitik sıfırlar yalnızca tekillik tespitiyle tanınır.

lambda_max=0 temasları -85.173420989 ve 94.746406171 derecededir.
Bunlar k=1 denklemi üzerinden fzero ile iyileştirilmiş tekil olmayan,
ancak SIKI contraction sağlamayan izole temas noktalarıdır; aralık sınırı değildir.

## Case seçimi ve özdeğerler

Case II, eski nominal çalışma noktasına en yakın kararlı tarama noktasıdır.
Case III, güç maksimumunun 30 derece ötesine en yakın pozitif güçlü kararsız
tarama noktasıdır. Bu yalnızca temsilci seçim kuralıdır; fiziksel parametreler
ayarlanmaz. Her biri için güç referansı denge denkleminden bulunur.

| Case | gamma0 (derece) | Pstar (pu) | eig(A) | eig(M), iki katlı | eig(Sym(F)) |
|---|---:|---:|---|---:|---|
| I | Mevcut değil | - | - | - | Sertifika yok |
| II | 9.555454969 | 1.000000000 | -41.29982645; -34.48824169 | 14120.5031934 | -750.58137010; +674.79330196 |
| III | 124.800000000 | 5.699033857 | -84.27186185; +8.48379371 | 1531770.8605581 | -397.86667377; +322.07860563 |

Case II yerel olarak kararlı fakat bu metrikle sertifikasızdır.
Case III pozitif bir A özdeğeri nedeniyle yerel olarak kararsızdır.
Her case'in iki M özdeğeri kendi içinde eşittir: M=beta0²*I. Grafiklerde
üst üste görünmeleri doğrudur. Farklı case'ler arasında değerler farklıdır.
Mutlak metrik büyüklükleri fiziksel kararlılık marjı değildir.

## Makaledeki Şekil 4 düzeni

`07_case_trajectories_figure4.png/.fig` iki satır, üç sütun içerir:
üstte faz düzlemleri, altta dengeye metrik uzaklık L(t).
Case I bulunmadığından onun iki panelinde matematiksel gerekçe gösterilir.
Case II ve III için gerçek nonlinear yörüngeler hesaplanır.

Denge yörüngesi sabittir. Pertürbe yörünge, dengeye M metriğinde
L(0)=0.2 uzaklıktan başlar. Case II'de başlangıç yönü Sym(F)'nin en büyük
özdeğer yönünden, Case III'te A'nın kararsız özvektöründen seçilir.
Bunlar başlangıç pertürbasyonlarıdır; denge hız farkı sıfır olsa da pertürbe
başlangıç hız farkı sıfır olmak zorunda değildir. Tam başlangıç durumları
`results.trajectories{2/3}.initialState` içindedir.

L(t)=||Theta_d*(x(t)-x_d)||, denge noktasında kurulan SABİT metrikle ölçülür.
Açı sarılmadan (unwrapped) kullanılır. Hız rad/s, zaman saniyedir.
Sertifika olmadığı için herhangi bir kanıtlanmış üstel azalma eğrisi çizilmez.

- Case II: 0.579908 s; L(0)=0.2, maksimum L yaklaşık 2.013,
  son L yaklaşık 5.998e-8. Geçici büyümenin ardından aynı dengeye döner.
- Case III: 2.35744 s; L(0)=0.2, maksimum L yaklaşık 42114.9,
  son L yaklaşık 6479.72. Başlangıçtaki kararsız dengeden uzaklaşır;
  sonlu zaman grafiği sınırsız büyüme iddiası değildir. Yörünge başka bir
  senkron çalışma noktasına yaklaşabilir.

Her iki simülasyon daha sıkı tolerans ve yarı maksimum adımla tekrar kontrol
edildi. Göreli farklar Case II için 2.62e-12, Case III için 2.18e-10.
Tüm taramada denge ve dönüşüm artıkları sıfır; doğrudan/analitik Sym(F)
maksimum farkı yaklaşık 1.61e-13. Denge dışındaki beta zaman türevleri de
bağımsız yönlü sonlu farkla kontrol edilir; esas hesaplar analitiktir.

## Dosyalar

- `baseline_gfm_contraction.m`: aynı, düzeltilmiş ana betik.
- `baseline_gfm_results.mat`: parametreler, tarama, sertifika teşhisi, case'ler,
  özdeğerler, dönüşüm doğrulamaları ve tüm yörünge verileri.
- `baseline_gfm_diagnostics.txt`: parametreler, denklemler, matrisler ve artıklar.
- `gfm_sg_case_comparison.csv`, `baseline_gfm_candidates.csv`: güncel case tablosu.
- `01_power_angle.png/.fig`: güç-açı bağıntısı ve gerçek çalışma noktaları.
- `02_contraction_interval.png/.fig`: boş contraction bölgesi ve sıfır temasları.
- `03_symF_eigenvalues.png/.fig`: iki Sym(F) özdeğeri.
- `04_gamma_omega_map.png/.fig`: bu sabit denge metriğiyle sertifikalı şerit yok.
- `05_case_eigenvalue_comparison.png/.fig`: üç panelli A, M, Sym(F) grafiği.
- `06_sg_operating_and_state_regions.png/.fig`: eski sınır grafiğinin düzeltilmiş hali;
  artık geçersiz olan eski dar aralıkları göstermez.
- `07_case_trajectories_figure4.png/.fig`: Şekil 4 mantığında altı panel.

Kaynak: Pourazad vd., “Large-signal strength assessment of grid connecting
assets using contraction analysis”, Bölüm 2.2, 2.3, 3.2 ve Şekil 3-4.
Kaynağın SG eğrileri GFM sonucuymuş gibi kopyalanmamıştır; doğru GFM
zaman türevleri nedeniyle oluşmayan kategori açıkça bildirilmiştir.
