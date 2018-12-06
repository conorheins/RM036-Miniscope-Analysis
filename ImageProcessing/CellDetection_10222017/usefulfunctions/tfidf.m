function tfidf_norm = tfidf(matrix)

% creates Term-Frequency Inverse-Document-Frequency transform on
% time-series data (T x N, where T is 'term' length and N is 'document'
% length


% checks for empty 'terms' and eliminates them (since TF term can't divide
% by 0)

if ~isempty(sum(matrix,2)==0)
    matrix = matrix(sum(matrix,2)>0,:);
end

tfidf_norm = zeros(size(matrix));
TF = bsxfun(@rdivide,matrix,sum(matrix,2));
IDF = log(bsxfun(@ldivide,sum(matrix),length(matrix)));
tfidf_norm = bsxfun(@times,TF,IDF);
tfidf_norm(isnan(tfidf_norm))=0;

end
