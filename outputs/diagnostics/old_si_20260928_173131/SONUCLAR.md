# Akım sınırlayıcısız VSM-GFM: dirençsiz SMIB

Aynı MATLAB dosyası yerinde güncellendi. Devre: E — jXeq — U1; E=U1=1, Xeq=0.3 pu. Pe=(E U1/Xeq) sin(gamma). Önceki dirençli ağ kullanılmıyor.

GFMVSG1 VSM parametreleri: H=0.1 s, K=20, Kd=40. İdeal PLL ve sabit frekans referansı 1 pu alınır. f0=50 Hz, Mv=2H/omega0=0.000636619772368, Dv=(K+Kd)/omega0=0.19098593171, Dv/Mv=300. Her gamma0 için Pstar=Pe(gamma0), Delta_omega0=0.

Nonlinear denklem: gamma_dot=w; w_dot=[(Pstar-Pe)/nu-Dv*w]/Mv, nu=1+w/omega0. 1/nu korunur. Normal hız bölgesi nu>0.5 dışında hesap durur; bir sınırlandırıcı eklenmez. Gerçek Simulink PLL, iç kontrol, filtre ve elektrik ağı dinamikleri bu iki durumlu modelde temsil edilmez.

## Yöntemlerin ayrımı

Ana yöntem: makaledeki beta_dot_paper=-C omega0 cos(gamma0) katsayısıyla her denge için ayrı Theta ve M tasarlanır. Bu katsayı gerçek VSM beta zaman türevi olarak sunulmaz. Theta_ref seçildikten sonra o dengeye ait analizde sabittir. Tam doğrudan F hesabında Theta_dot_ref=0 bu nedenle elde edilir; Eq15c/2 ile uyuşur.

Method B: beta=-C sin(gamma)/nu; beta_dot=-C cos(gamma) w/nu+C sin(gamma) f2/(omega0 nu^2). İkinci türev Hessian ve tam Jacobian ile hesaplanır. Gerçek dengede ikisi de sıfır olur ve bu yöntemin farklı metriği contraction sertifikası vermez.

Ayrıca gamma_dot=omega0 varsayımsal hareketi boyunca Theta türevi de hesaplanır. Bu hareket gerçek VSM denge yörüngesi değildir ve Eq15c ile uyuşmaz; ayrı tanı eğrisi olarak gösterilir. Böylece eşitlik yalnızca aynı sabit metrikle yapılan karşılaştırmada iddia edilir.

## Açılar ve aralıklar

SG kontrolü eski kod gibi gamma_rep=0, gamma0 ayrı olacak şekilde PASS verdi. GFM ana taraması yeni direktif uyarınca gerçek denge açısında gamma=gamma0 yapılır. gamma_rep=0 sonucu her case için ayrıca saklanır; ana sınıflandırma yerine kullanılmaz.

Ana paper-metric ve aynı metrikle doğrudan tanımın ortak gamma0 aralıkları:

- 47.4749771988° < gamma0 < 47.5586153232°; sınır güçleri 2.45660739, 2.45989367 pu.
- 89.9989592722° < gamma0 < 89.9996018686°; sınır güçleri 3.33333333, 3.33333333 pu.

Bunlar farklı dengelerde ayrı metriklerle elde edilen çalışma noktası aralıklarıdır; tek ortak metrikle global sertifika iddiası değildir.

| Case | gamma0 (deg) | Pstar (pu) | lambda max paper | lambda max direct |
|---|---:|---:|---:|---:|
| Case I | 47.516796261 | 2.458251181 | -12.291029772 | -12.291029772 |
| Case II | 17.457603124 | 1.000000000 | 351466.846135192 | 351466.846135192 |
| Case III | 105.000249983 | 3.219748990 | 15343.852960485 | 15343.852960485 |

Pstar<=1 pu sektörü gamma0<=17.45760312° ile sınırlı. Bulunan contracting dengeler 1 pu üstündedir; matematiksel sonuç olarak tutulmuş ve grafiklerde ayrılmıştır.

## Yörüngeler ve testler

Üç case için gerçek nonlinear model simüle edilir. Case I üstel sınırı yalnızca denge özdeğerinden alınmaz: tam nonlinear Jacobian için açı-hız dikdörtgeninde iki değişken giriş sınırlandırılır, köşe özdeğerleriyle contraction doğrulanır ve bu dikdörtgen içine yerleştirilen metrik topunun değişmezliği kullanılır. Önceki enerji sertifikası kullanılmaz. Case II/III için üstel azalma sınırı çizilmez.

Analitik A/B sonlu farkla, mutlak/sapma hız koordinatları eşdeğerlikle doğrulandı. Tüm case denge kalıntıları ve daha sıkı ODE toleransları kontrol edildi. MATLAB verileri baseline_gfm_results.mat, case tablosu gfm_sg_case_comparison.csv, ayrıntılar baseline_gfm_diagnostics.txt içindedir.
