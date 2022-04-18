from Cluster import Cluster, PartitionSet


class QualityMeasuerer:
    def calculate_quality(self, cluster: Cluster, partition_count: int) -> float:
        if len(cluster) == 0:
            return 0.0
        
        mean = cluster.mean_member_simularity(partition_count)
        sum_value = sum([ pow(cluster.member_simularity(item, partition_count) - mean, 2) for item in cluster])
        
        return sum_value / len(cluster)
    
    def calculate_speculative_quality(self, initial_quality: float, include_item: object, \
        cluster: Cluster, gamma: PartitionSet) -> float:
        if include_item in cluster:
            return initial_quality

        new_total_participation = len(gamma)+1
        new_mean = (cluster.sum_membership()+1) / (len(cluster)+1)
        sum_value = sum([pow((1 / new_total_participation) - new_mean, 2)] +\
            [pow(sim - new_mean, 2) for sim in cluster.calc_all_membersimularity(new_total_participation).values() ])

        return sum_value / (len(cluster)+1)
def target_bin_3_4th_count_estimator(gamma: PartitionSet) -> int:
    partition_ln = [len(partition) for partition in gamma]
    average = sum(partition_ln) / len(partition_ln)
    third = (max(partition_ln) - average ) / 2
    return int(average + third)