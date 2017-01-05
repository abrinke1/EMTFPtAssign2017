
Float_t GetMedian(const TH1D* hist);

void WeightByResScore(TH1D& hist, const Float_t med_ratio);

Float_t GetResScore(const TH1D* hist, const Float_t med_ratio);

Float_t GetResScoreErr(const TH1D* hist, const Float_t med_ratio);
